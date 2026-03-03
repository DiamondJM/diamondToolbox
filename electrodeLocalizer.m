classdef electrodeLocalizer < handle
% electrodeLocalizer  Localizes implanted electrodes from pre-op MRI and
% post-op CT, producing a leads.csv and cortical surface files that
% sourceLocalizer requires.
%
% USAGE
%   Constructed automatically by sourceLocalizer. Do not call directly
%   unless you need to re-run a specific stage.
%
%   sl = sourceLocalizer(subj, rootFolder);
%   sl = sourceLocalizer(subj, rootFolder, chanNames);
%   sl = sourceLocalizer(subj, rootFolder, chanNames, ...
%                        'forceNewElectrodeLocalizer', true);
%
% PIPELINE STAGES
%   1. setupDirectories  — create zloc/ folder structure under tal/
%   2. getInputFiles     — file dialogs for MRI + CT nifti; copy to zloc/
%   3. runSurface        — FreeSurfer recon-all + pial-outer-smoothed envelope
%   4. runSuma           — AFNI/SUMA standard ld141 mesh → gifti files
%   5. coregisterCT      — AFNI rigid-body CT→MR registration
%   6. detectElectrodes  — AFNI 3dclust on CT → cluster centroids in MR space
%   7. namingGUI         — interactive 3-D figure; user names each cluster,
%                          labels it depth or subdural, or marks as artifact
%   8. projectElectrodes — subdural contacts snapped to nearest pial vertex;
%                          depth contacts kept at CT-MR coordinates
%   9. writeLeads        — write tal/leads.csv (chanName, x, y, z)
%
% DEPENDENCIES
%   FreeSurfer  — recon-all surface reconstruction
%   AFNI/SUMA   — @SUMA_Make_Spec_FS standard mesh (ld141, 198 812 vertices)
%   AFNI        — align.sh CT→MR coregistration; 3dclust electrode detection
%   gifti       — bundled in Utilities/eeg_toolbox/localize/zLocalize/packages/
%
% OUTPUTS
%   tal/leads.csv                         — chanName, x, y, z (RAS mm)
%   tal/zloc/freesurfer/<subj>/SUMA/      — gifti surface files
%
% See also: sourceLocalizer, create_surf, suma

    properties
        subj
        rootFolder
        chanNames       % cell array of channel name strings (may be empty)

        locDirs         % struct of directory paths (from localizer_create_directories)

        % Detected electrode clusters (set by detectElectrodes)
        clusters        % [N x 3] centroid coordinates in MR mm space

        % Named leads (set by namingGUI + projectElectrodes)
        leads           % table with columns: chanName, x, y, z, type
    end

    properties (Constant)
        CT_HU_THRESHOLD          = 2500;  % HU threshold for electrode detection (Hounsfield units)
                                         %   1500: skull forms one huge connected blob (293K voxels)
                                         %         that swallows top-of-head contacts via guide bolts
                                         %   2000: skull separates but contact blooms touch/merge on
                                         %         close-spaced oblique shanks → fuses entire shank
                                         %   2500: good balance; micro-contact oblique shanks still
                                         %         merge (sub-voxel gap, CT resolution limit) but
                                         %         all other contacts resolve individually
        AFNI_CLUSTER_RMM         = 1;     % 3dclust: max neighbor distance in mm (rmm)
        AFNI_CLUSTER_VMUL        = 3;     % 3dclust: minimum cluster volume in µl (vmul)
        AFNI_CLUSTER_MAX_UL      = 1000;  % max cluster volume in µl; skull blob ~15,000 µl,
                                         %   guide bolts ~1,400 µl; contacts 5–50 µl
        SURFACE_PROXIMITY_MAX_MM = 40;    % max distance (mm) from pial surface to keep a cluster
        AFNI_MORPH_EROSION_ITER  = 3;    % morphological erosion iterations before clustering
                                         %   (each iter ≈ 1 voxel; 3 iters ≈ 1.25 mm at 0.415 mm CT)
                                         %   separates adjacent-contact blooms without destroying contacts
        AFNI_CLUSTER_VMUL_ERODED = 2;    % min cluster volume (µl) after erosion (contacts shrink ~10×)
        AFNI_BIN_DEFAULT         = fullfile(char(java.lang.System.getProperty('user.home')), 'abin');
    end

    properties (Hidden = true)
        fsBin    % FreeSurfer binary directory
        afniBin  % AFNI binary directory
    end

    methods

        % -----------------------------------------------------------------
        %% Constructor
        % -----------------------------------------------------------------

        function self = electrodeLocalizer(subj, rootFolder, chanNames, varargin)
            % electrodeLocalizer  Check completion and run pipeline if needed.
            %
            % Inputs:
            %   subj       - subject identifier string
            %   rootFolder - root data directory
            %   chanNames  - cell array of channel names (may be empty)
            %
            % Optional name-value:
            %   forceNew        - re-run even if complete (default false)
            %   freesurfer_bin  - path to FreeSurfer bin directory
            %   afni_bin        - path to AFNI bin directory

            p = inputParser;
            addParameter(p, 'forceNew',       false);
            addParameter(p, 'freesurfer_bin', electrodeLocalizer.defaultFsBin());
            addParameter(p, 'afni_bin',       electrodeLocalizer.AFNI_BIN_DEFAULT);
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;
            self.fsBin   = p.Results.freesurfer_bin;
            self.afniBin = p.Results.afni_bin;

            self.subj       = subj;
            self.rootFolder = rootFolder;
            self.chanNames  = chanNames;

            self.setupDirectories();

            if ~forceNew && self.isComplete()
                fprintf('[electrodeLocalizer] Localization already complete for %s. Skipping.\n', subj);
                return;
            end

            self.run('forceNew', forceNew);
        end

        % -----------------------------------------------------------------
        %% Completion check
        % -----------------------------------------------------------------

        function tf = isComplete(self)
            % Returns true if all 5 required localization files exist:
            %   tal/leads.csv, lh.pial.gii, rh.pial.gii,
            %   lh.pial-outer-smoothed.gii, rh.pial-outer-smoothed.gii

            leadsFile = fullfile(self.rootFolder, self.subj, 'tal', 'leads.csv');
            sumaDir   = fullfile(self.locDirs.fs_subj, 'SUMA');
            files = { ...
                leadsFile; ...
                fullfile(sumaDir, 'lh.pial.gii'); ...
                fullfile(sumaDir, 'rh.pial.gii'); ...
                fullfile(sumaDir, 'lh.pial-outer-smoothed.gii'); ...
                fullfile(sumaDir, 'rh.pial-outer-smoothed.gii'); ...
            };
            tf = all(cellfun(@(f) exist(f,'file')==2, files));
        end

        % -----------------------------------------------------------------
        %% Top-level pipeline orchestrator
        % -----------------------------------------------------------------

        function run(self, varargin)
            % Run the full localization pipeline from imaging to leads.csv.
            %
            % Optional name-value:
            %   forceNew - re-run all stages (default false); individual
            %              stages always check their own outputs first.

            p = inputParser;
            addParameter(p, 'forceNew', false);
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;

            fprintf('\n+----------------------------------------------------------+\n');
            fprintf('|            electrodeLocalizer — %s\n', self.subj);
            fprintf('+----------------------------------------------------------+\n');

            % Show setup dialog until complete, cancelled, or user picks Create.
            while ~self.isComplete()
                dlg = self.localizationSetupDialog();
                switch dlg.action
                    case 'cancel'
                        error('electrodeLocalizer:cancelled', ...
                            '[electrodeLocalizer] Setup cancelled by user.');
                    case 'import'
                        % Dialog copied whatever the user provided; loop back
                        % so they can see updated status and import the rest,
                        % or exit automatically if everything is now in place.
                        continue;
                    case 'create'
                        break;   % proceed to full pipeline below
                end
            end

            % If import made everything complete, nothing left to do.
            if self.isComplete()
                return;
            end

            self.checkPrerequisites('errorIfMissing', true);
            self.getInputFiles();
            self.runSurface();           % never force — recon-all takes hours
            self.runSuma();              % never force — SUMA takes minutes
            self.coregisterCT('forceNew', forceNew);
            self.detectElectrodes('forceNew', forceNew);
            if isempty(self.chanNames)
                self.chanNames = sourceLocalizer.loadChanNamesFromFile();
            end
            self.namingGUI();
            self.projectElectrodes();
            self.writeLeads();
        end

        % -----------------------------------------------------------------
        %% Import existing localization files
        % -----------------------------------------------------------------

        function importExistingLocalization(self)
            % Copy existing leads.csv and/or gifti surface files into the
            % expected folder locations, replacing the need to run part or
            % all of the localization pipeline.
            %
            % Each file is handled independently: if it is already present
            % at the expected location it is reported and skipped; if it is
            % missing a file dialog is opened. Files that the user dismisses
            % are skipped silently so a partial import is always valid.
            %
            % Can be called standalone at any time:
            %   sl.electrodeLocalizer.importExistingLocalization()
            %
            % Validates:
            %   leads.csv  — must contain columns chanName, x, y, z.
            %                Warns if any chanNames entries are absent.
            %   *.gii      — must be loadable as a gifti object.

            fprintf('\n[electrodeLocalizer] Import existing localization files for %s.\n', self.subj);
            fprintf('Dismiss any dialog to skip that file.\n\n');

            talDir  = fullfile(self.rootFolder, self.subj, 'tal');
            sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');
            if ~exist(sumaDir, 'dir'), mkdir(sumaDir); end

            % --- leads.csv ---
            leadsFile = fullfile(talDir, 'leads.csv');
            if exist(leadsFile, 'file') == 2
                fprintf('[import] leads.csv already present: %s\n', leadsFile);
            else
                fprintf('Select leads.csv (columns required: chanName, x, y, z)...\n');
                [f, d] = uigetfile('*.csv', 'Select leads.csv');
                if ~isequal(f, 0)
                    src = fullfile(d, f);
                    electrodeLocalizer.validateLeadsCSV(src, self.chanNames);
                    if ~exist(talDir, 'dir'), mkdir(talDir); end
                    copyfile(src, leadsFile);
                    fprintf('[import] leads.csv copied to %s\n', leadsFile);
                else
                    fprintf('[import] leads.csv skipped.\n');
                end
            end

            % --- Gifti surfaces ---
            surfaceFiles = { ...
                'lh.pial-outer-smoothed.gii', 'left  pial-outer-smoothed (required for projection)'; ...
                'rh.pial-outer-smoothed.gii', 'right pial-outer-smoothed (required for projection)'; ...
                'lh.pial.gii',                'left  pial (used for electrode naming display)';       ...
                'rh.pial.gii',                'right pial (used for electrode naming display)';       ...
            };

            for i = 1:size(surfaceFiles, 1)
                fname = surfaceFiles{i, 1};
                label = surfaceFiles{i, 2};
                dest  = fullfile(sumaDir, fname);

                if exist(dest, 'file') == 2
                    fprintf('[import] %s already present.\n', fname);
                    continue;
                end

                fprintf('Select %s — %s...\n', fname, label);
                [f, d] = uigetfile('*.gii', sprintf('Select %s', fname));
                if isequal(f, 0)
                    fprintf('[import] %s skipped.\n', fname);
                    continue;
                end

                src = fullfile(d, f);
                electrodeLocalizer.validateGifti(src);
                copyfile(src, dest);
                fprintf('[import] %s copied to %s\n', fname, sumaDir);
            end

            fprintf('\n[electrodeLocalizer] Import complete.\n');
            if self.isComplete()
                fprintf('All required files are now in place.\n\n');
            else
                fprintf('Note: some required files are still missing.\n');
                fprintf('Call checkPrerequisites() for a full status report.\n\n');
            end
        end

        % -----------------------------------------------------------------
        %% Prerequisites check
        % -----------------------------------------------------------------

        function prereqs = checkPrerequisites(self, varargin)
            % Check all external dependencies and report their status.
            % Only flags an item as REQUIRED if the corresponding pipeline
            % stage has not yet produced its output files.
            %
            % Usage:
            %   sl.electrodeLocalizer.checkPrerequisites()      % report only
            %   sl.electrodeLocalizer.checkPrerequisites('errorIfMissing', true)
            %
            % Returns a struct with logical fields:
            %   .imageProcessingToolbox
            %   .gifti
            %   .freesurfer
            %   .afni

            p = inputParser;
            addParameter(p, 'errorIfMissing', false);
            parse(p, varargin{:});
            errorIfMissing = p.Results.errorIfMissing;

            % Determine which stages still need to run
            surfDir    = fullfile(self.locDirs.fs_subj, 'surf');
            sumaDir    = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhPial     = fullfile(surfDir, 'lh.pial-outer-smoothed');
            rhPial     = fullfile(surfDir, 'rh.pial-outer-smoothed');
            lhGii      = fullfile(sumaDir, 'lh.pial-outer-smoothed.gii');
            rhGii      = fullfile(sumaDir, 'rh.pial-outer-smoothed.gii');
            xfmFile    = fullfile(self.locDirs.ct_1_xfm, 'transform.mat');
            clustFile  = fullfile(self.locDirs.ct_1_xfm, 'clusters_mr.mat');

            needSurface = ~(exist(lhPial,'file')==2 && exist(rhPial,'file')==2);
            needSuma    = ~(exist(lhGii,'file')==2  && exist(rhGii,'file')==2);
            needCoreg   = exist(xfmFile, 'file') ~= 2;
            needDetect  = exist(clustFile,'file') ~= 2;

            % --- Check each prerequisite ---

            % gifti toolbox
            prereqs.gifti = exist('gifti','file') == 2;

            % FreeSurfer — check that recon-all exists at the expected bin path.
            % Executing it directly from MATLAB fails because recon-all is a
            % shell script that requires FREESURFER_HOME to be set, which GUI
            % MATLAB does not inherit from the user's shell profile.
            fsCmd = fullfile(self.fsBin, 'recon-all');
            prereqs.freesurfer = exist(fsCmd, 'file') == 2;

            % AFNI/SUMA — check that @SUMA_Make_Spec_FS exists at the expected bin path.
            afniCmd = fullfile(self.afniBin, '@SUMA_Make_Spec_FS');
            prereqs.afni = exist(afniCmd, 'file') == 2;

            % --- Determine which missing items are actually required ---
            missing = {};

            if ~prereqs.gifti
                missing{end+1} = sprintf( ...
                    '  [REQUIRED] gifti toolbox\n%s\n%s', ...
                    '             Needed for: surface display and projection (stages 7-8)', ...
                    '             Fix: run addpath(genpath(pwd)) from the diamondToolbox root.');
            end

            if ~prereqs.freesurfer && needSurface
                missing{end+1} = sprintf( ...
                    '  [REQUIRED] FreeSurfer\n%s\n%s\n%s', ...
                    '             Needed for: surface reconstruction (stage 3)', ...
                    sprintf('             Expected bin: %s', self.fsBin), ...
                    '             Fix: install FreeSurfer or pass ''freesurfer_bin'' to sourceLocalizer.');
            end

            if ~prereqs.afni && needSuma
                missing{end+1} = sprintf( ...
                    '  [REQUIRED] AFNI/SUMA\n%s\n%s\n%s', ...
                    '             Needed for: standard surface mesh (stage 4)', ...
                    sprintf('             Expected bin: %s', self.afniBin), ...
                    '             Fix: install AFNI or pass ''afni_bin'' to sourceLocalizer.');
            end

            % --- Print report ---
            fprintf('\n+----------------------------------------------------------+\n');
            fprintf('|       electrodeLocalizer — Prerequisites (%s)\n', self.subj);
            fprintf('+----------------------------------------------------------+\n');
            electrodeLocalizer.printPrereqLine('gifti toolbox',             prereqs.gifti,                 true);
            electrodeLocalizer.printPrereqLine(sprintf('FreeSurfer  (%s)', self.fsBin),  prereqs.freesurfer, needSurface);
            electrodeLocalizer.printPrereqLine(sprintf('AFNI/SUMA   (%s)', self.afniBin), prereqs.afni,      needSuma);
            fprintf('+----------------------------------------------------------+\n');

            if isempty(missing)
                fprintf('| All required prerequisites satisfied.\n');
            else
                fprintf('| %d prerequisite(s) missing:\n\n', numel(missing));
                for i = 1:numel(missing)
                    fprintf('%s\n\n', missing{i});
                end
            end
            fprintf('+----------------------------------------------------------+\n\n');

            if errorIfMissing && ~isempty(missing)
                error('[electrodeLocalizer] Missing prerequisites. See report above.');
            end
        end

        % -----------------------------------------------------------------
        %% Stage 1 — directory setup
        % -----------------------------------------------------------------

        function setupDirectories(self)
            self.locDirs = localizer_create_directories(self.subj, self.rootFolder);
        end

        % -----------------------------------------------------------------
        %% Stage 2 — acquire input image files
        % -----------------------------------------------------------------

        function getInputFiles(self)
            % Prompt for pre-op MRI and post-op CT imaging files if not
            % already present in the zloc folder structure.
            %
            % Accepted formats: .nii, .nii.gz, .mgz
            %   .nii     — copied directly.
            %   .nii.gz  — decompressed via MATLAB gunzip.
            %   .mgz     — converted via FreeSurfer mri_convert (uses fsBin).

            filter = {'*.nii;*.nii.gz;*.mgz', 'Imaging files (*.nii, *.nii.gz, *.mgz)'};

            mrDest  = fullfile(self.locDirs.mr_pre, 'mr_pre.nii');
            ctDest  = fullfile(self.locDirs.ct_1,   'ct_implant.nii');

            % MRI
            if exist(mrDest, 'file') ~= 2
                while true
                    choice = dlgNonModal( ...
                        {'Select the pre-operative T1 MPRAGE MRI for this subject.', '', ...
                         'Accepted formats:  .nii  |  .nii.gz  |  .mgz'}, ...
                        'Pre-op MRI', 'Browse...', 'Cancel');
                    if ~strcmp(choice, 'Browse...')
                        error('[electrodeLocalizer] MRI is required.');
                    end
                    [f, d] = uigetfile(filter, 'Select pre-op MRI');
                    if ~isequal(f, 0), break; end
                    % file picker cancelled — loop back to description dialog
                end
                self.convertToNii(fullfile(d, f), mrDest);
            else
                fprintf('[Stage 2] MRI already present: %s\n', mrDest);
            end

            % CT
            if exist(ctDest, 'file') ~= 2
                while true
                    choice = dlgNonModal( ...
                        {'Select the post-operative CT scan for this subject.', '', ...
                         'This should be the CT acquired after electrode implantation.', ...
                         'Accepted formats:  .nii  |  .nii.gz  |  .mgz'}, ...
                        'Post-op CT', 'Browse...', 'Cancel');
                    if ~strcmp(choice, 'Browse...')
                        error('[electrodeLocalizer] CT is required.');
                    end
                    [f, d] = uigetfile(filter, 'Select post-op CT');
                    if ~isequal(f, 0), break; end
                    % file picker cancelled — loop back to description dialog
                end
                self.convertToNii(fullfile(d, f), ctDest);
            else
                fprintf('[Stage 2] CT already present: %s\n', ctDest);
            end
        end

        % -----------------------------------------------------------------
        %% Stage 3 — FreeSurfer surface reconstruction
        % -----------------------------------------------------------------

        function runSurface(self, varargin)
            p = inputParser;
            addParameter(p, 'forceNew', false);
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;

            surfDir = fullfile(self.locDirs.fs_subj, 'surf');
            lhPial  = fullfile(surfDir, 'lh.pial-outer-smoothed');
            rhPial  = fullfile(surfDir, 'rh.pial-outer-smoothed');
            nativeDone = exist(lhPial, 'file') == 2 && exist(rhPial, 'file') == 2;

            % Also treat imported SUMA GIFTIs as sufficient — the user may
            % have surfaces from a prior FreeSurfer run and just imported them.
            sumaDir  = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhGii    = fullfile(sumaDir, 'lh.pial-outer-smoothed.gii');
            rhGii    = fullfile(sumaDir, 'rh.pial-outer-smoothed.gii');
            giftiDone = exist(lhGii, 'file') == 2 && exist(rhGii, 'file') == 2;

            if (nativeDone || giftiDone) && ~forceNew
                fprintf('[Stage 3] Surface outputs already present; skipping recon-all.\n');
                return;
            end

            fprintf('[Stage 3] Running FreeSurfer surface reconstruction...\n');

            % Pre-set FREESURFER (bare, no _HOME) so recon-all can find it
            % even when its own `setenv LANG C` (issued before sourcing
            % FreeSurferEnv.csh) prevents the source chain from setting it.
            % Also ensures FREESURFER_HOME is in the environment for any
            % FreeSurfer tool that needs it.
            fsHome = fileparts(self.fsBin);
            setenv('FREESURFER',      fsHome);
            setenv('FREESURFER_HOME', fsHome);

            mrNii = fullfile(self.locDirs.mr_pre, 'mr_pre.nii');
            create_surf(self.subj, mrNii, self.locDirs.fs, ...
                'freesurfer_bin', self.fsBin);
        end

        % -----------------------------------------------------------------
        %% Stage 4 — SUMA standard mesh
        % -----------------------------------------------------------------

        function runSuma(self, varargin)
            p = inputParser;
            addParameter(p, 'forceNew', false);
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;

            sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhGii   = fullfile(sumaDir, 'lh.pial-outer-smoothed.gii');
            rhGii   = fullfile(sumaDir, 'rh.pial-outer-smoothed.gii');
            done    = exist(lhGii, 'file') == 2 && exist(rhGii, 'file') == 2;

            if done && ~forceNew
                fprintf('[Stage 4] SUMA surfaces already exist; skipping.\n');
                return;
            end

            fprintf('[Stage 4] Running SUMA to generate standard ld141 mesh...\n');

            % @SUMA_Make_Spec_FS requires mri/orig.mgz and an mri/orig/
            % directory (even if empty). Subjects from public datasets often
            % omit these; T1.mgz is equivalent to orig.mgz.
            mriDir  = fullfile(self.locDirs.fs_subj, 'mri');
            origMgz = fullfile(mriDir, 'orig.mgz');
            origDir = fullfile(mriDir, 'orig');
            t1Mgz   = fullfile(mriDir, 'T1.mgz');
            if exist(origMgz, 'file') ~= 2 && exist(t1Mgz, 'file') == 2
                fprintf('[Stage 4] orig.mgz not found; copying T1.mgz as orig.mgz.\n');
                copyfile(t1Mgz, origMgz);
            end
            if exist(origDir, 'dir') ~= 7
                fprintf('[Stage 4] mri/orig/ directory not found; creating it.\n');
                mkdir(origDir);
            end

            suma(self.subj, self.locDirs.fs, ...
                'afni_bin',       self.afniBin, ...
                'freesurfer_bin', self.fsBin, ...
                'rerun',          forceNew);

            % @SUMA_Make_Spec_FS sometimes exits 0 even on failure.
            % Verify the expected output actually exists.
            if exist(lhGii, 'file') ~= 2 || exist(rhGii, 'file') ~= 2
                error('[Stage 4] SUMA ran but expected output not found:\n  %s\n  %s', ...
                    lhGii, rhGii);
            end
        end

        % -----------------------------------------------------------------
        %% Stage 5 — CT-to-MR coregistration via AFNI
        % -----------------------------------------------------------------

        function coregisterCT(self, varargin)
            % Rigid-body coregistration of post-op CT to pre-op MRI using
            % AFNI align.sh (LPC cost function, @Align_Centers initialisation).
            % Saves the path to the .aff12.1D combined transform to
            % zloc/CT_1/transform/transform.mat.

            p = inputParser;
            addParameter(p, 'forceNew', false);
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;

            xfmFile = fullfile(self.locDirs.ct_1_xfm, 'transform.mat');
            if exist(xfmFile, 'file') == 2 && ~forceNew
                fprintf('[Stage 5] CT-MR transform already exists; skipping.\n');
                return;
            end

            assert(exist(fullfile(self.afniBin, 'align_epi_anat.py'), 'file') == 2, ...
                '[electrodeLocalizer] AFNI not found at %s. Pass ''afni_bin'' to sourceLocalizer.', self.afniBin);

            fprintf('[Stage 5] Coregistering CT to MR via AFNI align.sh...\n');

            mrNii = fullfile(self.locDirs.mr_pre, 'mr_pre.nii');
            ctNii = fullfile(self.locDirs.ct_1,   'ct_implant.nii');
            assert(exist(mrNii,'file')==2, '[electrodeLocalizer] MRI not found: %s', mrNii);
            assert(exist(ctNii,'file')==2, '[electrodeLocalizer] CT not found: %s',  ctNii);

            alignScript = fullfile(fileparts(mfilename('fullpath')), ...
                'Utilities', 'eeg_toolbox', 'localize', 'zLocalize', ...
                'shell_scripts', 'align.sh');
            assert(exist(alignScript,'file')==2, ...
                '[electrodeLocalizer] align.sh not found at: %s', alignScript);

            workDir = self.locDirs.ct_1_xfm;
            % Copy inputs to work dir so align.sh can find them by stem.
            % Pre-orient the MR NIfTI to RAI before align.sh runs.
            % align.sh calls "3dWarp -deoblique" which has no -overwrite flag
            % and, when applied to an RSP-oriented NIfTI with negative deltas,
            % shifts the physical coordinate origin relative to orig.mgz/SUMA
            % space (causing a systematic ~20-35 mm electrode displacement).
            % Pre-reorienting to RAI with 3dresample makes 3dWarp a no-op,
            % keeping the coordinate frame consistent with the SUMA surfaces.
            mrWork = fullfile(workDir, 'mr_pre.nii');
            ctWork = fullfile(workDir, 'ct_implant.nii');
            if exist(mrWork, 'file') ~= 2
                fprintf('[Stage 5] Resampling MR to RAI orientation for AFNI registration...\n');
                cmd = sprintf('"%s" -orient RAI -inset "%s" -prefix "%s" -overwrite', ...
                    fullfile(self.afniBin, '3dresample'), mrNii, mrWork);
                [cst, ctx] = unix(cmd);
                if cst ~= 0
                    fprintf('%s\n', ctx);
                    fprintf('[Stage 5] 3dresample failed; copying mr_pre.nii directly.\n');
                    copyfile(mrNii, mrWork);
                end
            end
            if exist(ctWork,'file') ~= 2, copyfile(ctNii, ctWork); end

            setenv('PATH', [getenv('PATH') ':' self.afniBin]);

            % Skip align.sh if the combined transform already exists.
            % align_epi_anat.py can hang on its 3dNotes history step even after
            % the actual alignment finishes; skipping avoids repeated hangs.
            hits = dir(fullfile(workDir, 'full_*.aff12.1D'));
            if isempty(hits)
                % Delete stale align.sh intermediates so a retry runs clean.
                % align.sh tools (3dWarp, 3dAutomask, @Align_Centers, etc.)
                % don't accept -overwrite and will silently skip or error on
                % existing files, causing align_epi_anat.py to run on stale data.
                stale = [dir(fullfile(workDir, 'mr_pre_do*')); ...
                         dir(fullfile(workDir, 'ct_implant_shft*')); ...
                         dir(fullfile(workDir, 'ct_implant_sh2*')); ...
                         dir(fullfile(workDir, '__tt_*')); ...
                         dir(fullfile(workDir, '*_mat.aff12.1D'))];
                for si = 1:numel(stale)
                    delete(fullfile(stale(si).folder, stale(si).name));
                end

                % Run align.sh in the background so the known 3dNotes hang
                % (AFNI 26.x) does not block MATLAB indefinitely.
                % align_epi_anat.py writes *_mat.aff12.1D BEFORE calling 3dNotes,
                % so we poll for that file instead of waiting for shell exit.
                logFile = fullfile(workDir, 'align_log.txt');
                unix(sprintf('bash "%s" mr_pre ct_implant "%s" > "%s" 2>&1 &', ...
                    alignScript, workDir, logFile));

                fprintf('[Stage 5] Running align_epi_anat.py (polling every 15 s) .');
                t0_align = tic; alignTimeout = 3600;
                xfmAf   = []; xfmShft = fullfile(workDir, 'ct_implant_shft.1D');
                while toc(t0_align) < alignTimeout
                    xfmAf = dir(fullfile(workDir, '*_XFMTO_lpc_*_mat.aff12.1D'));
                    if ~isempty(xfmAf) && exist(xfmShft, 'file') == 2
                        break;
                    end
                    pause(15);
                    fprintf('.');
                end
                fprintf(' (%.0f min)\n', toc(t0_align)/60);

                % Kill any hanging 3dNotes / align_epi_anat children.
                unix('pkill -f "3dNotes" 2>/dev/null; pkill -f "align_epi_anat" 2>/dev/null');

                % Build the combined transform from the per-step files.
                hits = dir(fullfile(workDir, 'full_*.aff12.1D'));
                if isempty(hits) && ~isempty(xfmAf) && exist(xfmShft,'file') == 2
                    xfmAfPath   = fullfile(workDir, xfmAf(1).name);
                    xfmFullPath = fullfile(workDir, ['full_' xfmAf(1).name]);
                    fprintf('[Stage 5] Running cat_matvec to combine transforms...\n');
                    catCmd = sprintf('cat_matvec -ONELINE "%s" "%s" > "%s"', ...
                        xfmAfPath, xfmShft, xfmFullPath);
                    [cstatus, ctxt] = unix(catCmd);
                    if cstatus ~= 0, fprintf('%s\n', ctxt); end
                    hits = dir(fullfile(workDir, 'full_*.aff12.1D'));
                end

                assert(~isempty(hits), ...
                    '[Stage 5] Alignment did not produce a .aff12.1D transform.\n  workDir: %s', workDir);
            else
                fprintf('[Stage 5] Existing AFNI transform found; skipping align.sh.\n');
            end

            aff1D = fullfile(workDir, hits(1).name);
            fprintf('[Stage 5] Using AFNI transform: %s\n', hits(1).name);

            save(xfmFile, 'aff1D');
            fprintf('[Stage 5] Transform path saved to %s\n', xfmFile);
        end

        % -----------------------------------------------------------------
        %% Stage 6 — electrode cluster detection via AFNI
        % -----------------------------------------------------------------

        function detectElectrodes(self, varargin)
            % Detect electrode contact positions on the co-registered CT BRIK.
            % Centroids are stored in self.clusters as [N x 3] MR-space mm.
            %
            % Pipeline:
            %   1. 3dclust on raw CT at CT_HU_THRESHOLD → clusters + clst.1D
            %      (skull forms one large separate cluster at 2000 HU — no masking needed)
            %   2. Volume filter: remove clusters > AFNI_CLUSTER_MAX_UL (removes skull/hardware)
            %   3. xfm_leads: transform CT centroids → MR space via AFNI
            %   4. Proximity filter: drop clusters > SURFACE_PROXIMITY_MAX_MM from pial

            p = inputParser;
            addParameter(p, 'forceNew', false);
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;

            clustFile = fullfile(self.locDirs.ct_1_xfm, 'clusters_mr.mat');
            if exist(clustFile, 'file') == 2 && ~forceNew
                fprintf('[Stage 6] Cluster file already exists; loading.\n');
                S = load(clustFile, 'clusters_mm');
                self.clusters = S.clusters_mm;
                return;
            end

            fprintf('[Stage 6] Detecting electrode clusters via AFNI 3dclust...\n');

            % Load the AFNI transform path saved by coregisterCT
            xfmFile = fullfile(self.locDirs.ct_1_xfm, 'transform.mat');
            assert(exist(xfmFile, 'file') == 2, ...
                '[electrodeLocalizer] Run coregisterCT before detectElectrodes.');
            S = load(xfmFile, 'aff1D');
            aff1D = S.aff1D;

            % Verify the .aff12.1D file still exists; fall back to glob search
            if exist(aff1D, 'file') ~= 2
                hits = dir(fullfile(self.locDirs.ct_1_xfm, 'full_*.aff12.1D'));
                assert(~isempty(hits), ...
                    '[Stage 6] AFNI transform (.aff12.1D) not found. Re-run coregisterCT.');
                aff1D = fullfile(self.locDirs.ct_1_xfm, hits(1).name);
                fprintf('[Stage 6] Using transform: %s\n', hits(1).name);
            end

            % CT BRIK created by align.sh (3dresample -orient RAI).
            % AFNI writes compressed BRIK.gz by default; check both forms.
            workDir = self.locDirs.ct_1_xfm;
            ctBrik  = fullfile(workDir, 'ct_implant+orig');
            assert(exist([ctBrik '.BRIK'], 'file') == 2 || ...
                   exist([ctBrik '.BRIK.gz'], 'file') == 2, ...
                '[Stage 6] CT BRIK not found: %s.BRIK(.gz)\n  Re-run coregisterCT.', ctBrik);

            % Ensure AFNI tools are on PATH
            setenv('PATH', [getenv('PATH') ':' self.afniBin]);

            % Clustering subdirectory
            clustDir = fullfile(workDir, 'clustering');
            if ~exist(clustDir, 'dir'), mkdir(clustDir); end

            % Parameters
            rmm  = electrodeLocalizer.AFNI_CLUSTER_RMM;          % neighbor dist in mm
            vmul = electrodeLocalizer.AFNI_CLUSTER_VMUL_ERODED;  % min volume (µl) — small for high-HU blobs
            cv   = electrodeLocalizer.CT_HU_THRESHOLD;            % HU threshold

            outNii  = sprintf('3dclusters_r%d_is%d_thr%d.nii', rmm, vmul, cv);
            outFull = fullfile(clustDir, outNii);
            clst1D  = fullfile(clustDir, 'clst.1D');

            fprintf('[Stage 6] Parameters: rmm=%dmm, vmul=%dµl, threshold=%dHU\n', rmm, vmul, cv);

            % Run 3dclust directly on the raw CT (no skull masking).
            % At 2000 HU the skull forms one enormous cluster (~15,000 µl) that is
            % removed by AFNI_CLUSTER_MAX_UL.  Skull masking was tried but the
            % skull-stripped brain mask doesn't cover the superior brain apex, so
            % it zeroes out top-of-head contacts — causing them to go undetected.
            cmd1 = sprintf( ...
                '3dclust -savemask "%s" -overwrite -1Dformat -1clip %d %d %d "%s" > "%s"', ...
                outFull, cv, rmm, vmul, ctBrik, clst1D);
            [status, txt] = unix(sprintf('cd "%s" && %s', clustDir, cmd1));
            if status ~= 0
                fprintf('%s\n', txt);
                error('[Stage 6] 3dclust failed (exit %d).', status);
            end

            % Parse clst.1D: centroids (CT RAI mm) + cluster volumes.
            % Filter by volume: 3dclust's vmul already removes too-small blobs;
            % AFNI_CLUSTER_MAX_UL removes skull/bone artifacts (typically > 2000 µl).
            [centroids_ct, volumes] = electrodeLocalizer.parseClst1D(clst1D);
            keep = volumes <= electrodeLocalizer.AFNI_CLUSTER_MAX_UL;
            n_total   = size(centroids_ct, 1);
            n_removed = sum(~keep);
            centroids_ct = centroids_ct(keep, :);
            N = size(centroids_ct, 1);
            fprintf('[Stage 6] Found %d clusters; removed %d artifact(s) > %d µl; keeping %d.\n', ...
                n_total, n_removed, electrodeLocalizer.AFNI_CLUSTER_MAX_UL, N);
            assert(N > 0, ...
                '[Stage 6] No clusters found. Check CT_HU_THRESHOLD (%d HU) or CT file.', cv);

            % Step 4: transform CT RAI mm → MR space via xfm_leads
            % (uses @ElectroGrid + ConvertSurface from AFNI, same pipeline as ALICE3)
            chanNums  = arrayfun(@(i) sprintf('e%03d',i), (1:N)', 'UniformOutput', false);
            ct_coords = table(chanNums, ...
                centroids_ct(:,1), centroids_ct(:,2), centroids_ct(:,3), ...
                'VariableNames', {'chanName','x','y','z'});

            mr_coords = xfm_leads(ct_coords, aff1D, ctBrik, 'RAI', clustDir);
            clusters_mm = mr_coords{:, {'x','y','z'}};

            % Spatial proximity filter: drop clusters too far from the pial
            % surface.  Hardware anchors and dense bone fragments that survive
            % the volume filter often land > 40 mm from the nearest cortical
            % vertex.  Depth electrode contacts in hippocampus/amygdala are
            % typically 20-35 mm from cortex, so a 40 mm threshold keeps them.
            sumaDir  = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhPialF  = fullfile(sumaDir, 'lh.pial.gii');
            rhPialF  = fullfile(sumaDir, 'rh.pial.gii');
            if exist(lhPialF,'file') == 2 && exist(rhPialF,'file') == 2
                pialV = [gifti(lhPialF).vertices; gifti(rhPialF).vertices];
                maxD  = electrodeLocalizer.SURFACE_PROXIMITY_MAX_MM;
                prox_keep = false(size(clusters_mm,1), 1);
                for ci = 1:size(clusters_mm,1)
                    d = clusters_mm(ci,:) - pialV;          % Nv × 3
                    prox_keep(ci) = sqrt(min(sum(d.^2,2))) <= maxD;
                end
                n_prox = sum(~prox_keep);
                clusters_mm = clusters_mm(prox_keep, :);
                fprintf('[Stage 6] Removed %d cluster(s) > %d mm from pial surface; %d remaining.\n', ...
                    n_prox, maxD, size(clusters_mm,1));
            end

            self.clusters = clusters_mm;
            save(clustFile, 'clusters_mm');
            fprintf('[Stage 6] %d cluster centroids saved to %s\n', size(clusters_mm,1), clustFile);
        end

        % -----------------------------------------------------------------
        %% Stage 7 — interactive naming GUI
        % -----------------------------------------------------------------

        function namingGUI(self)
            % Single-window electrode naming GUI.
            %
            % Left panel  : rotatable 3-D brain with all clusters plotted.
            %               Drag to rotate at any time — no dialogs block it.
            % Right panel : channel list, free-text field, type selector,
            %               Confirm / Artifact / Back / Finish buttons.
            %
            % Results stored in self.leads (chanName, x, y, z, type).

            assert(~isempty(self.clusters), ...
                '[electrodeLocalizer] Run detectElectrodes before namingGUI.');

            N = size(self.clusters, 1);
            fprintf('[Stage 7] Launching electrode naming GUI (%d clusters)...\n', N);

            % Reorder clusters by greedy nearest-neighbor, starting from the
            % most-anterior point (max y in RAS).  This keeps contacts on the
            % same shank consecutive and groups nearby shanks together, making
            % the naming task much faster than the arbitrary volume-sort order.
            pts = self.clusters;
            order   = zeros(1, N);
            visited = false(1, N);
            [~, seed] = max(pts(:, 2));   % most anterior cluster first
            order(1) = seed;  visited(seed) = true;
            for ii = 2:N
                prev = pts(order(ii-1), :);
                d = sum((pts - prev).^2, 2);
                d(visited) = Inf;
                [~, nxt] = min(d);
                order(ii) = nxt;  visited(nxt) = true;
            end
            localClusters = pts(order, :);   % N×3, greedy-NN order

            % Sort chanNames: alphabetical by shank tag, then numerical by
            % contact number within each tag (e.g. LACN 1…8, LAMY 1…7, …).
            sortedChanNames = self.chanNames;
            if ~isempty(sortedChanNames)
                toks = regexp(sortedChanNames, '^(.*)\s+(\d+)$', 'tokens', 'once');
                valid_tok = ~cellfun(@isempty, toks);
                keys = sortedChanNames;  % fallback: sort as-is
                keys(valid_tok) = cellfun(@(t) sprintf('%s_%05d', t{1}, str2double(t{2})), ...
                    toks(valid_tok), 'UniformOutput', false);
                [~, si] = sort(keys);
                sortedChanNames = sortedChanNames(si);
            end

            sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhFile  = fullfile(sumaDir, 'lh.pial.gii');
            rhFile  = fullfile(sumaDir, 'rh.pial.gii');
            assert(exist(lhFile,'file')==2 && exist(rhFile,'file')==2, ...
                '[electrodeLocalizer] SUMA pial .gii not found; re-run Stage 4.');
            lhSurf = gifti(lhFile);
            rhSurf = gifti(rhFile);

            % ---- mutable state shared by callbacks (closures) ----
            k          = 1;
            nameOut    = cell(N, 1);
            typeOut    = cell(N, 1);
            artifact   = false(N, 1);

            % ---- figure ----
            fig = figure('Name', sprintf('Electrode Naming — %s', self.subj), ...
                'NumberTitle', 'off', 'Color', [0.12 0.12 0.12], ...
                'Position', [50 50 1400 820], ...
                'CloseRequestFcn', @cbClose);

            % left: 3-D brain axes (63% width, full height)
            ax = axes('Parent', fig, 'Position', [0.01 0.01 0.62 0.97], ...
                'Color', 'k', 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
            hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
            view(ax,3); camlight(ax,'headlight'); material(ax,'dull');

            patch(ax, 'Faces', lhSurf.faces, 'Vertices', lhSurf.vertices, ...
                'FaceColor', [0.75 0.70 0.65], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.4, 'PickableParts', 'none', 'HitTest', 'off');
            patch(ax, 'Faces', rhSurf.faces, 'Vertices', rhSurf.vertices, ...
                'FaceColor', [0.75 0.70 0.65], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.4, 'PickableParts', 'none', 'HitTest', 'off');

            dotColors = repmat([0 0.9 0.9], N, 1);   % cyan — visible against brain surface
            hDots = scatter3(ax, localClusters(:,1), localClusters(:,2), ...
                localClusters(:,3), 50, dotColors, 'filled', 'HitTest', 'off');
            hHi  = scatter3(ax, NaN, NaN, NaN, 140, [1 0.8 0], 'filled', ...
                'MarkerEdgeColor', 'w', 'LineWidth', 1.5, 'HitTest', 'off');

            rotate3d(ax, 'on');   % always-on rotation — drag to rotate

            % ---- right panel controls ----
            bg  = [0.18 0.18 0.18];
            fg  = [0.92 0.92 0.92];
            x0  = 0.645;   w = 0.345;

            hInfo = uicontrol(fig, 'Style', 'text', ...
                'Units','normalized','Position',[x0 0.88 w 0.10], ...
                'BackgroundColor',bg,'ForegroundColor',fg, ...
                'FontSize',13,'FontWeight','bold','HorizontalAlignment','center');

            hXYZ = uicontrol(fig, 'Style', 'text', ...
                'Units','normalized','Position',[x0 0.83 w 0.05], ...
                'BackgroundColor',bg,'ForegroundColor',[0.5 0.8 1.0], ...
                'FontSize',10,'HorizontalAlignment','center');

            uicontrol(fig,'Style','text','String','Channel name:', ...
                'Units','normalized','Position',[x0 0.79 w 0.04], ...
                'BackgroundColor',bg,'ForegroundColor',fg, ...
                'FontSize',10,'HorizontalAlignment','left');

            hList = uicontrol(fig, 'Style', 'listbox', ...
                'Units','normalized','Position',[x0 0.52 w 0.27], ...
                'BackgroundColor',[0.22 0.22 0.22],'ForegroundColor',fg, ...
                'FontSize',11,'String',{});

            uicontrol(fig,'Style','text','String','Or type name:', ...
                'Units','normalized','Position',[x0 0.47 w 0.04], ...
                'BackgroundColor',bg,'ForegroundColor',fg, ...
                'FontSize',10,'HorizontalAlignment','left');

            hEdit = uicontrol(fig, 'Style', 'edit', ...
                'Units','normalized','Position',[x0 0.42 w 0.05], ...
                'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',fg, ...
                'FontSize',11,'HorizontalAlignment','left','String','');

            uicontrol(fig,'Style','text','String','Electrode type:', ...
                'Units','normalized','Position',[x0 0.37 w 0.04], ...
                'BackgroundColor',bg,'ForegroundColor',fg, ...
                'FontSize',10,'HorizontalAlignment','left');

            hType = uicontrol(fig, 'Style', 'popupmenu', ...
                'Units','normalized','Position',[x0 0.32 w 0.05], ...
                'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',fg, ...
                'FontSize',11,'String',{'depth','subdural'},'Value',1);

            uicontrol(fig,'Style','pushbutton', ...
                'String','Confirm & Next  >', ...
                'Units','normalized','Position',[x0 0.24 w 0.07], ...
                'BackgroundColor',[0.18 0.48 0.18],'ForegroundColor','w', ...
                'FontSize',12,'FontWeight','bold','Callback',@cbConfirm);

            uicontrol(fig,'Style','pushbutton', ...
                'String','Mark as Artifact', ...
                'Units','normalized','Position',[x0 0.16 w 0.07], ...
                'BackgroundColor',[0.48 0.18 0.18],'ForegroundColor','w', ...
                'FontSize',12,'Callback',@cbArtifact);

            uicontrol(fig,'Style','pushbutton','String','< Back', ...
                'Units','normalized','Position',[x0 0.08 w*0.44 0.07], ...
                'BackgroundColor',[0.30 0.30 0.30],'ForegroundColor',fg, ...
                'FontSize',11,'Callback',@cbBack);

            uicontrol(fig,'Style','pushbutton','String','Finish', ...
                'Units','normalized','Position',[x0+w*0.56 0.08 w*0.44 0.07], ...
                'BackgroundColor',[0.20 0.38 0.58],'ForegroundColor','w', ...
                'FontSize',11,'FontWeight','bold','Callback',@cbFinish);

            uicontrol(fig,'Style','text', ...
                'String','Drag on brain to rotate at any time', ...
                'Units','normalized','Position',[x0 0.02 w 0.05], ...
                'BackgroundColor',bg,'ForegroundColor',[0.50 0.50 0.50], ...
                'FontSize',9,'HorizontalAlignment','center');

            refreshDisplay();
            waitfor(fig);   % non-blocking for controls; blocks here until closed

            % ---- build leads table after GUI closes ----
            valid    = ~artifact & ~cellfun(@isempty, nameOut);
            nameOut  = nameOut(valid);
            typeOut  = typeOut(valid);
            xyz      = localClusters(valid, :);
            self.leads = table(nameOut, xyz(:,1), xyz(:,2), xyz(:,3), typeOut, ...
                'VariableNames', {'chanName','x','y','z','type'});
            fprintf('[Stage 7] %d contacts accepted, %d artifacts.\n', ...
                sum(valid), sum(artifact));

            % ---- nested callbacks (share k, nameOut, etc. by closure) ----

            function refreshDisplay()
                if ~ishandle(fig), return; end
                set(hInfo, 'String', sprintf('Cluster  %d  /  %d', k, N));
                set(hXYZ,  'String', sprintf('x=%.1f   y=%.1f   z=%.1f', ...
                    localClusters(k,:)));
                set(hHi, 'XData', localClusters(k,1), ...
                         'YData', localClusters(k,2), ...
                         'ZData', localClusters(k,3));

                % List box: remaining (unassigned) channel names, sorted
                used = nameOut(~cellfun(@isempty, nameOut));
                if ~isempty(sortedChanNames)
                    remaining = sortedChanNames(~ismember(sortedChanNames, used));
                    set(hList, 'String', remaining(:), ...
                        'Value', max(1, min(get(hList,'Value'), numel(remaining))));
                else
                    set(hList, 'String', {});
                end

                % Restore prior entry when navigating back
                set(hEdit, 'String', '');
                if ~isempty(nameOut{k})
                    set(hEdit, 'String', nameOut{k});
                end
                typeStrs = get(hType, 'String');
                if ~isempty(typeOut{k})
                    idx = find(strcmp(typeStrs, typeOut{k}), 1);
                    if ~isempty(idx), set(hType, 'Value', idx); end
                end
                drawnow;
            end

            function name = resolvedName()
                name = strtrim(get(hEdit, 'String'));
                if ~isempty(name), return; end
                listStr = get(hList, 'String');
                if ~isempty(listStr)
                    name = listStr{max(1, get(hList,'Value'))};
                end
            end

            function cbConfirm(~,~)
                name = resolvedName();
                if isempty(name)
                    msgbox('Enter a channel name or use Mark as Artifact.','','warn');
                    return;
                end
                typeStrs    = get(hType,'String');
                nameOut{k}  = name;
                typeOut{k}  = typeStrs{get(hType,'Value')};
                artifact(k) = false;
                c = get(hDots,'CData');
                if strcmp(typeOut{k},'depth'), col=[0.3 0.6 1.0]; else, col=[0.3 1.0 0.5]; end
                c(k,:) = col;
                set(hDots,'CData',c);
                if k < N, k = k+1; end
                refreshDisplay();
            end

            function cbArtifact(~,~)
                artifact(k) = true;
                nameOut{k}  = '';
                typeOut{k}  = '';
                c = get(hDots,'CData'); c(k,:)=[0.6 0.2 0.2]; set(hDots,'CData',c);
                if k < N, k = k+1; end
                refreshDisplay();
            end

            function cbBack(~,~)
                if k > 1
                    k = k - 1;
                    % Undo: clear the naming for the cluster we just returned to,
                    % reset its dot color to cyan, and let refreshDisplay put the
                    % name back into the available list automatically.
                    nameOut{k}  = '';
                    typeOut{k}  = '';
                    artifact(k) = false;
                    c = get(hDots, 'CData');
                    c(k, :) = [0 0.9 0.9];
                    set(hDots, 'CData', c);
                end
                refreshDisplay();
            end

            function cbFinish(~,~)
                delete(fig);
            end

            function cbClose(~,~)
                delete(fig);
            end

        end

        % -----------------------------------------------------------------
        %% Stage 8 — project subdural contacts to pial surface
        % -----------------------------------------------------------------

        function projectElectrodes(self)
            % Depth contacts: keep CT-MR coordinates unchanged.
            % Subdural contacts: snap to nearest vertex on the ipsilateral
            %   pial-outer-smoothed surface (brain-shift correction without
            %   inter-electrode spacing constraints).

            assert(~isempty(self.leads), ...
                '[electrodeLocalizer] Run namingGUI before projectElectrodes.');

            sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhFile  = fullfile(sumaDir, 'lh.pial-outer-smoothed.gii');
            rhFile  = fullfile(sumaDir, 'rh.pial-outer-smoothed.gii');
            lhSurf  = gifti(lhFile);
            rhSurf  = gifti(rhFile);

            xyz = self.leads{:, {'x','y','z'}};

            for i = 1:height(self.leads)
                if ~strcmp(self.leads.type{i}, 'subdural')
                    continue;   % depth: no projection
                end

                % Hemisphere by x-coordinate sign (RAS: x<0 = left)
                if xyz(i, 1) < 0
                    verts = lhSurf.vertices;
                else
                    verts = rhSurf.vertices;
                end

                % Nearest vertex on the envelope surface
                dists = sum((verts - xyz(i,:)).^2, 2);
                [~, idx] = min(dists);
                xyz(i, :) = verts(idx, :);
            end

            self.leads.x = xyz(:, 1);
            self.leads.y = xyz(:, 2);
            self.leads.z = xyz(:, 3);

            fprintf('[Stage 8] Subdural contacts projected to pial surface.\n');
        end

        % -----------------------------------------------------------------
        %% Stage 9 — write leads.csv
        % -----------------------------------------------------------------

        function writeLeads(self)
            % Write tal/leads.csv with columns chanName, x, y, z.
            % Coordinates are rounded to 2 decimal places.

            assert(~isempty(self.leads), ...
                '[electrodeLocalizer] Run projectElectrodes before writeLeads.');

            talDir    = fullfile(self.rootFolder, self.subj, 'tal');
            if ~exist(talDir, 'dir'), mkdir(talDir); end
            leadsFile = fullfile(talDir, 'leads.csv');

            out      = self.leads(:, {'chanName','x','y','z'});
            out.x    = round(out.x, 2);
            out.y    = round(out.y, 2);
            out.z    = round(out.z, 2);

            writetable(out, leadsFile);
            fprintf('[Stage 9] leads.csv written to %s\n', leadsFile);
        end

    end % methods

    methods (Access = private)

        function convertToNii(self, srcFile, destFile)
            % Convert an imaging file to uncompressed NIfTI at destFile.
            %
            % Supported inputs:
            %   .nii     — copied directly.
            %   .nii.gz  — decompressed with MATLAB gunzip.
            %   .mgz     — converted via FreeSurfer's MRIread/MRIwrite MATLAB tools.

            if endsWith(lower(srcFile), '.nii.gz')
                tmpDir   = tempname;
                mkdir(tmpDir);
                result   = gunzip(srcFile, tmpDir);
                movefile(result{1}, destFile);
                rmdir(tmpDir, 's');
            else
                [~, ~, ext] = fileparts(srcFile);
                switch lower(ext)
                    case '.nii'
                        copyfile(srcFile, destFile);
                    case '.mgz'
                        % Use FreeSurfer's bundled MATLAB tools to read .mgz
                        % and write .nii — avoids mri_convert binary entirely
                        % (which may be ARM64-native while MATLAB runs x86_64
                        % under Rosetta on Apple Silicon).
                        fsMatlabDir = fullfile(fileparts(self.fsBin), 'matlab');
                        if exist(fsMatlabDir, 'dir') ~= 7
                            error('[electrodeLocalizer] FreeSurfer MATLAB tools not found at: %s', fsMatlabDir);
                        end
                        addpath(fsMatlabDir);
                        mri = MRIread(srcFile);
                        err = MRIwrite(mri, destFile);
                        if err ~= 0
                            error('[electrodeLocalizer] MRIwrite failed writing: %s', destFile);
                        end
                    otherwise
                        error('[electrodeLocalizer] Unsupported format: %s', ext);
                end
            end
        end

    end % methods (Access = private)

    methods (Static, Access = private)

        function p = defaultFsBin()
            % Locate the FreeSurfer bin directory by checking, in order:
            %   1. FREESURFER_HOME environment variable
            %   2. Version subdirectories under /Applications/freesurfer/
            %      (handles installs like /Applications/freesurfer/8.1.0/)
            %   3. /Applications/freesurfer/bin directly
            % Whichever location first contains recon-all is returned.

            base = '/Applications/freesurfer';

            % 1. Env var (works when MATLAB is launched from a shell)
            fsHome = getenv('FREESURFER_HOME');
            if ~isempty(fsHome) && ...
                    exist(fullfile(fsHome, 'bin', 'recon-all'), 'file') == 2
                p = fullfile(fsHome, 'bin');
                return;
            end

            % 2. Version subdirectories — sort descending so latest wins
            if exist(base, 'dir') == 7
                d = dir(base);
                d = d([d.isdir] & ~strncmp({d.name}, '.', 1));
                names = fliplr(sort({d.name}));
                for i = 1:numel(names)
                    candidate = fullfile(base, names{i}, 'bin');
                    if exist(fullfile(candidate, 'recon-all'), 'file') == 2
                        p = candidate;
                        return;
                    end
                end

                % 3. Direct install (no version subdir)
                if exist(fullfile(base, 'bin', 'recon-all'), 'file') == 2
                    p = fullfile(base, 'bin');
                    return;
                end
            end

            % Fallback — let checkPrerequisites report the miss clearly
            p = fullfile(base, 'bin');
        end

        function validateLeadsCSV(filepath, chanNames)
            % Validate that a CSV file is a usable leads table.
            % Errors if required columns are absent.
            % Warns if chanNames entries are missing from the file.

            try
                T = readtable(filepath);
            catch e
                error('[electrodeLocalizer] Could not read %s: %s', filepath, e.message);
            end

            required = {'chanName','x','y','z'};
            missing  = required(~ismember(required, T.Properties.VariableNames));
            if ~isempty(missing)
                error('[electrodeLocalizer] leads.csv is missing required columns: %s', ...
                    strjoin(missing, ', '));
            end

            if ~isempty(chanNames)
                absent = setdiff(chanNames, T.chanName);
                if ~isempty(absent)
                    warning('[electrodeLocalizer] %d channel(s) in chanNames not found in leads.csv:\n  %s', ...
                        numel(absent), strjoin(absent, ', '));
                end
            end
        end

        function validateGifti(filepath)
            % Validate that a file is a loadable gifti object.
            try
                g = gifti(filepath);  %#ok<NASGU>
            catch e
                error('[electrodeLocalizer] Could not load gifti file %s: %s', ...
                    filepath, e.message);
            end
        end

        function [centroids, volumes] = parseClst1D(clst1DPath)
            % Parse AFNI 3dclust -1Dformat output and return cluster centroids + volumes.
            %
            % Returns:
            %   centroids - [N x 3] matrix of centroid coordinates in CT RAI mm
            %               (CM_RL, CM_AP, CM_IS — columns 2, 3, 4 of each data row)
            %   volumes   - [N x 1] vector of cluster volumes in µl (column 1)

            fid = fopen(clst1DPath, 'r');
            assert(fid > 0, '[electrodeLocalizer] Cannot open %s', clst1DPath);
            lines = {};
            while ~feof(fid)
                ln = strtrim(fgetl(fid));
                if ischar(ln) && ~isempty(ln) && ln(1) ~= '#'
                    lines{end+1} = ln; %#ok<AGROW>
                end
            end
            fclose(fid);

            if isempty(lines)
                centroids = zeros(0, 3);
                volumes   = zeros(0, 1);
                return;
            end

            % Each data row: Volume  CM_RL  CM_AP  CM_IS  minRL ...
            % Column 1 = Volume (µl); Columns 2, 3, 4 = RAI mm centroids.
            centroids = zeros(numel(lines), 3);
            volumes   = zeros(numel(lines), 1);
            for i = 1:numel(lines)
                vals = sscanf(lines{i}, '%f');
                volumes(i)     = vals(1);
                centroids(i,:) = vals(2:4)';
            end
        end

    end % methods (Static, Access = private) — utility functions

    methods (Access = private)

        function dlg = localizationSetupDialog(self)
            % Dark-themed modal dialog shown when required localization files
            % are missing.  A listbox shows all 5 files with ✓/✗ status;
            % selecting a row updates the description panel.  Create runs the
            % full pipeline; Import opens sequential file dialogs for each
            % missing file then copies them into place.
            %
            % Returns struct with .action: 'create' | 'import' | 'cancel'

            % ---- file list -----------------------------------------------
            talDir  = fullfile(self.rootFolder, self.subj, 'tal');
            sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');

            names = { ...
                'leads.csv'; ...
                'lh.pial.gii'; ...
                'rh.pial.gii'; ...
                'lh.pial-outer-smoothed.gii'; ...
                'rh.pial-outer-smoothed.gii'; ...
            };
            dests = { ...
                fullfile(talDir,  'leads.csv'); ...
                fullfile(sumaDir, 'lh.pial.gii'); ...
                fullfile(sumaDir, 'rh.pial.gii'); ...
                fullfile(sumaDir, 'lh.pial-outer-smoothed.gii'); ...
                fullfile(sumaDir, 'rh.pial-outer-smoothed.gii'); ...
            };
            descs = { ...
                sprintf(['CSV table of electrode coordinates.\n' ...
                    'Required columns: chanName, x, y, z.\n' ...
                    'Each row is one implanted contact in FreeSurfer RAS (mm).\n' ...
                    'Destination: %s'], dests{1}); ...
                sprintf(['Left hemisphere pial surface (GIFTI).\n' ...
                    'Generated by AFNI/SUMA (@SUMA_Make_Spec_FS).\n' ...
                    'Used for 3-D display in the electrode naming GUI.\n' ...
                    'Destination: %s'], dests{2}); ...
                sprintf(['Right hemisphere pial surface (GIFTI).\n' ...
                    'Generated by AFNI/SUMA (@SUMA_Make_Spec_FS).\n' ...
                    'Used for 3-D display in the electrode naming GUI.\n' ...
                    'Destination: %s'], dests{3}); ...
                sprintf(['Left pial-outer-smoothed surface (GIFTI).\n' ...
                    'Convex-hull envelope of the left pial surface.\n' ...
                    'Required for snapping subdural contacts to the cortex.\n' ...
                    'Destination: %s'], dests{4}); ...
                sprintf(['Right pial-outer-smoothed surface (GIFTI).\n' ...
                    'Convex-hull envelope of the right pial surface.\n' ...
                    'Required for snapping subdural contacts to the cortex.\n' ...
                    'Destination: %s'], dests{5}); ...
            };
            N       = numel(names);
            present = cellfun(@(f) exist(f,'file')==2, dests);

            % Listbox strings: status icon + filename
            listStrs = cell(N, 1);
            for i = 1:N
                if present(i)
                    listStrs{i} = ['[OK]  ' names{i}];
                else
                    listStrs{i} = ['[  ]  ' names{i}];
                end
            end

            % ---- theme ---------------------------------------------------
            BG   = [0.13 0.13 0.13];
            FG   = [0.92 0.92 0.92];
            DIM  = [0.52 0.52 0.52];
            BTN  = [0.25 0.25 0.25];
            EDBG = [0.20 0.20 0.20];
            LBBG = [0.18 0.18 0.18];

            % ---- layout --------------------------------------------------
            W      = 560;
            PAD    = 12;
            HDR_H  = 54;
            LST_H  = N * 28;
            DESC_H = 100;
            BTN_H  = 50;
            SEP    = 8;
            totalH = HDR_H + LST_H + SEP + DESC_H + SEP + BTN_H;

            % y positions (bottom = 0)
            btnY  = 8;
            descY = btnY  + BTN_H  + SEP;
            lstY  = descY + DESC_H + SEP;
            hdrY  = lstY  + LST_H;

            % ---- figure --------------------------------------------------
            ss  = get(0, 'ScreenSize');
            fig = figure( ...
                'Name',           sprintf('electrodeLocalizer — %s', self.subj), ...
                'MenuBar',        'none', ...
                'ToolBar',        'none', ...
                'NumberTitle',    'off', ...
                'Color',          BG, ...
                'Resize',         'off', ...
                'Position',       [round((ss(3)-W)/2) round((ss(4)-totalH)/2) W totalH], ...
                'CloseRequestFcn', @cbCancel);

            % ---- header --------------------------------------------------
            uicontrol(fig, 'Style','text', ...
                'String', sprintf('Localization files not found for  %s', self.subj), ...
                'ForegroundColor', FG, 'BackgroundColor', BG, ...
                'FontSize', 12, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'Position', [PAD hdrY+28 W-PAD*2 22]);
            uicontrol(fig, 'Style','text', ...
                'String', 'Select a row to see details.  Create runs the full pipeline.  Import copies existing files.', ...
                'ForegroundColor', DIM, 'BackgroundColor', BG, ...
                'FontSize', 9, 'HorizontalAlignment', 'left', ...
                'Position', [PAD hdrY+10 W-PAD*2 16]);

            % ---- description box (created first so hList callback can ref it) --
            uicontrol(fig, 'Style','text', 'String', 'File details', ...
                'ForegroundColor', DIM, 'BackgroundColor', BG, ...
                'FontSize', 9, 'HorizontalAlignment', 'left', ...
                'Position', [PAD descY+DESC_H-2 80 14]);
            hDesc = uicontrol(fig, 'Style','edit', ...
                'String', descs{1}, ...
                'ForegroundColor', FG, 'BackgroundColor', EDBG, ...
                'FontSize', 10, 'HorizontalAlignment', 'left', ...
                'Max', 5, 'Min', 0, 'Enable', 'inactive', ...
                'Position', [PAD descY W-PAD*2 DESC_H-16]);

            % ---- file listbox --------------------------------------------
            hList = uicontrol(fig, 'Style','listbox', ...
                'String', listStrs, ...
                'Value', 1, ...
                'BackgroundColor', LBBG, ...
                'ForegroundColor', FG, ...
                'FontSize', 11, ...
                'FontName', 'Courier', ...
                'Position', [PAD lstY W-PAD*2 LST_H], ...
                'Callback', @(src,~) set(hDesc, 'String', descs{get(src,'Value')})); %#ok<NASGU>

            % ---- action buttons ------------------------------------------
            bW = 130; bH = 34;
            uicontrol(fig, 'Style','pushbutton', ...
                'String', 'Create', ...
                'ForegroundColor', FG, 'BackgroundColor', [0.18 0.42 0.18], ...
                'FontSize', 11, 'FontWeight', 'bold', ...
                'Position', [PAD btnY bW bH], ...
                'Callback', @cbCreate);
            uicontrol(fig, 'Style','pushbutton', ...
                'String', 'Import', ...
                'ForegroundColor', FG, 'BackgroundColor', BTN, ...
                'FontSize', 11, ...
                'Position', [PAD + bW + PAD btnY bW bH], ...
                'Callback', @cbImport);
            uicontrol(fig, 'Style','pushbutton', ...
                'String', 'Cancel', ...
                'ForegroundColor', DIM, 'BackgroundColor', BTN, ...
                'FontSize', 11, ...
                'Position', [W - PAD - bW btnY bW bH], ...
                'Callback', @cbCancel);

            % ---- block until closed --------------------------------------
            setappdata(0, 'eloc_dlg_result', struct('action','cancel'));
            uiwait(fig);
            dlg = getappdata(0, 'eloc_dlg_result');
            if isappdata(0, 'eloc_dlg_result')
                rmappdata(0, 'eloc_dlg_result');
            end

            % ==== nested callbacks ========================================

            function cbCreate(~,~)
                setappdata(0, 'eloc_dlg_result', struct('action','create'));
                delete(fig);
            end

            function cbImport(~,~)
                % Open a file dialog for each missing file in turn.
                if ~exist(talDir,  'dir'), mkdir(talDir);  end
                if ~exist(sumaDir, 'dir'), mkdir(sumaDir); end
                for k = 1:N
                    if present(k), continue; end
                    if strcmp(names{k}, 'leads.csv')
                        filt = {'*.csv', 'CSV file (*.csv)'};
                    else
                        filt = {'*.gii', 'GIFTI surface (*.gii)'};
                    end
                    [f, d] = uigetfile(filt, sprintf('Select %s  (%d of %d missing)', ...
                        names{k}, sum(~present), N));
                    figure(fig);   % restore focus after uigetfile
                    if isequal(f, 0)
                        return;   % cancelled — leave dialog open
                    end
                    src = fullfile(d, f);
                    try
                        if strcmp(names{k}, 'leads.csv')
                            electrodeLocalizer.validateLeadsCSV(src, self.chanNames);
                        else
                            electrodeLocalizer.validateGifti(src);
                        end
                    catch e
                        warndlg(e.message, sprintf('Validation failed — %s', names{k}));
                        continue;
                    end
                    copyfile(src, dests{k});
                    fprintf('[import] %s → %s\n', src, dests{k});
                end
                setappdata(0, 'eloc_dlg_result', struct('action','import'));
                delete(fig);
            end

            function cbCancel(~,~)
                setappdata(0, 'eloc_dlg_result', struct('action','cancel'));
                if ishandle(fig), delete(fig); end
            end
        end

    end % methods (Access = private) — localizationSetupDialog

    methods (Static, Access = private)

        function out = ternary(cond, a, b)
            % Inline ternary: returns a if cond is true, else b.
            if cond, out = a; else, out = b; end
        end

        function printPrereqLine(label, found, needed)
            % Print one formatted line of the prerequisites report.
            %   found  — logical, whether the tool was detected
            %   needed — logical, whether this stage still needs to run
            if found
                status = 'OK     ';
            elseif needed
                status = 'MISSING';
            else
                status = 'OK     ';   % not found but not needed — treat as fine
            end
            if found
                tag = '';
            elseif needed
                tag = '  <-- required';
            else
                tag = '  (not needed — stage already complete)';
            end
            fprintf('| [%s]  %s%s\n', status, label, tag);
        end

    end % static private methods

end % classdef
