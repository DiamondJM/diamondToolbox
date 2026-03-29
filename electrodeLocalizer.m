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
        CT_HU_THRESHOLD          = 1800;  % HU threshold for electrode detection (Hounsfield units)
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
            self.fsBin      = p.Results.freesurfer_bin;
            self.afniBin    = p.Results.afni_bin;

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
            %   manual   - if true, replace auto-detect + naming GUI with
            %              the CT slicer manual localization GUI (default false).

            p = inputParser;
            addParameter(p, 'forceNew', false);
            addParameter(p, 'manual',   true);
            parse(p, varargin{:});
            forceNew  = p.Results.forceNew;
            useManual = p.Results.manual;

            fprintf('\n+----------------------------------------------------------+\n');
            fprintf('|            electrodeLocalizer — %s\n', self.subj);
            fprintf('+----------------------------------------------------------+\n');

            % Short-circuit if everything is already in place.
            if ~forceNew && self.isComplete(), return; end

            % Show setup dialog.  Imports are handled inline; dialog only
            % closes when the user clicks Create or Cancel.
            dlg = self.localizationSetupDialog();
            if strcmp(dlg.action, 'cancel')
                error('electrodeLocalizer:cancelled', ...
                    '[electrodeLocalizer] Setup cancelled by user.');
            end

            % If imports made everything complete, nothing left to do.
            if self.isComplete(), return; end

            self.checkPrerequisites('errorIfMissing', true);
            self.getInputFiles();
            self.runSurface();           % never force — recon-all takes hours
            self.runSuma();              % never force — SUMA takes minutes
            self.coregisterCT('forceNew', forceNew);
            if isempty(self.chanNames)
                self.chanNames = sourceLocalizer.loadChanNamesFromFile();
            end
            if useManual
                self.manualLocalize();
            else
                self.detectElectrodes('forceNew', forceNew);
                self.namingGUI();
            end
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

            filter = {'*.*', 'All files (*.nii, *.nii.gz, *.mgz)'};

            mrDest  = fullfile(self.locDirs.mr_pre, 'mr_pre.nii');
            ctDest  = fullfile(self.locDirs.ct_1,   'ct_implant.nii');

            % MRI
            if exist(mrDest, 'file') ~= 2
                while true
                    choice = dlgNonModal( ...
                        {'Select the pre-operative T1 MPRAGE MRI for this subject.', ...
                         'Ideally 1 mm isotropic — thicker slices reduce surface reconstruction quality.', ...
                         '', 'Accepted formats:  .nii  |  .nii.gz  |  .mgz'}, ...
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

            fprintf('[Stage 3] Running FreeSurfer surface reconstruction (recon-all — expect 8–10+ hours)...\n');

            % Pre-set FREESURFER (bare, no _HOME) so recon-all can find it
            % even when its own `setenv LANG C` (issued before sourcing
            % FreeSurferEnv.csh) prevents the source chain from setting it.
            % Also ensures FREESURFER_HOME is in the environment for any
            % FreeSurfer tool that needs it.
            fsHome = fileparts(self.fsBin);
            setenv('FREESURFER',      fsHome);
            setenv('FREESURFER_HOME', fsHome);

            mrNii = fullfile(self.locDirs.mr_pre, 'mr_pre.nii');

            rawPial   = fullfile(self.locDirs.fs_subj, 'surf', 'lh.pial');
            reconDone = exist(rawPial, 'file') == 2;
            fsSubjDir = self.locDirs.fs_subj;
            subjHasFiles = exist(fsSubjDir, 'dir') == 7 && numel(dir(fsSubjDir)) > 2;

            % If a prior partial run exists, resume it with -no-isrunning
            % rather than deleting it.  Deleting would throw away hours of
            % completed recon-all stages.  recon-all -no-isrunning -all
            % skips already-finished steps and continues from where it left off.
            isRunningFile = fullfile(self.locDirs.fs_subj, 'scripts', 'IsRunning.lh+rh');
            if subjHasFiles && ~reconDone
                if exist(isRunningFile, 'file') == 2
                    fprintf('[Stage 3] Stale IsRunning lock detected; removing so recon-all can resume.\n');
                    delete(isRunningFile);
                end
                fprintf('[Stage 3] Partial FreeSurfer output detected; resuming recon-all from last completed stage.\n');
            end

            if subjHasFiles && ~reconDone
                % Resume: pass -no-isrunning flag so recon-all doesn't abort
                % due to the (now deleted) lock file, and skips done steps.
                setenv('SUBJECTS_DIR', self.locDirs.fs);
                resumeCmd = sprintf('source %s && recon-all -sd %s -subjid %s -all -no-isrunning', ...
                    fullfile(fileparts(self.fsBin), 'SetUpFreeSurfer.sh'), ...
                    self.locDirs.fs, self.subj);
                fprintf('[Stage 3] Running: %s\n', resumeCmd);
                [st, txt] = unix(['bash -c ''' strrep(resumeCmd,'''','''''') '''']);
                if st ~= 0
                    fprintf('%s', txt);
                    error('[Stage 3] recon-all resume failed.');
                end
                % After resume, run envelope step
                create_surf(self.subj, mrNii, self.locDirs.fs, ...
                    'freesurfer_bin', self.fsBin, 'envelope_only', true);
            else
                envelopeOnly = reconDone && ~forceNew;
                create_surf(self.subj, mrNii, self.locDirs.fs, ...
                    'freesurfer_bin', self.fsBin, 'envelope_only', envelopeOnly);
            end
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

            % Ensure FREESURFER_HOME is set — required by mri_convert and
            % mris_convert called internally by @SUMA_Make_Spec_FS.
            fsHome = fileparts(self.fsBin);
            setenv('FREESURFER',      fsHome);
            setenv('FREESURFER_HOME', fsHome);

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
                'rerun',          ~done || forceNew);

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

            ctNii = fullfile(self.locDirs.ct_1, 'ct_implant.nii');
            assert(exist(ctNii,'file')==2, '[electrodeLocalizer] CT not found: %s', ctNii);
            mrNii      = fullfile(self.locDirs.mr_pre, 'mr_pre.nii');
            mrWorkName = 'mr_pre.nii';
            assert(exist(mrNii,'file')==2, '[electrodeLocalizer] MRI not found: %s', mrNii);

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
            mrWork = fullfile(workDir, mrWorkName);
            ctWork = fullfile(workDir, 'ct_implant.nii');
            if exist(mrWork, 'file') ~= 2
                fprintf('[Stage 5] Resampling MR to RAI orientation for AFNI registration...\n');
                cmd = sprintf('"%s" -orient RAI -inset "%s" -prefix "%s" -overwrite', ...
                    fullfile(self.afniBin, '3dresample'), mrNii, mrWork);
                [cst, ctx] = unix(cmd);
                if cst ~= 0
                    fprintf('%s\n', ctx);
                    fprintf('[Stage 5] 3dresample failed; copying %s directly.\n', mrWorkName);
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
                mrStem = mrWorkName(1:end-4);  % strip .nii → 'mr_pre' or 'mr_post'
                stale = [dir(fullfile(workDir, [mrStem '_do*'])); ...
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
                unix(sprintf('bash "%s" %s ct_implant "%s" > "%s" 2>&1 &', ...
                    alignScript, mrStem, workDir, logFile));

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
        %% Stage 6+7 (manual) — CT slicer + 3-D review
        % -----------------------------------------------------------------

        function manualLocalize(self)
            % Replaces detectElectrodes (Stage 6) + namingGUI (Stage 7).
            % Presents a three-panel CT slice viewer; user clicks contacts,
            % names each one, then reviews placements in 3-D.
            % Results stored in self.leads (chanName, x, y, z, type).

            ctFile = fullfile(self.locDirs.ct_1_xfm, 'ct_implant.nii');
            assert(exist(ctFile,'file')==2, ...
                '[manualLocalize] Coregistered CT not found: %s\nRun coregisterCT first.', ctFile);

            % Load CT→MR transform (produced by coregisterCT) for coord xfm.
            xfmFilePath = fullfile(self.locDirs.ct_1_xfm, 'transform.mat');
            assert(exist(xfmFilePath,'file')==2, ...
                '[manualLocalize] CT-MR transform not found. Run coregisterCT first.');
            S_xfm   = load(xfmFilePath, 'aff1D');
            aff1D   = S_xfm.aff1D;
            workDir = self.locDirs.ct_1_xfm;
            ctBrik  = fullfile(workDir, 'ct_implant+orig');
            assert(exist([ctBrik '.BRIK'],'file')==2 || exist([ctBrik '.BRIK.gz'],'file')==2, ...
                '[manualLocalize] CT BRIK not found: %s.BRIK(.gz)\nRun coregisterCT first.', ctBrik);

            % ---- Load CT ----
            info = niftiinfo(ctFile);
            vol  = double(niftiread(info));          % [nx, ny, nz]
            Txfm = info.Transform.T;                 % 4×4 MATLAB affine (1-based vox → world)
            [nx, ny, nz] = size(vol);

            % vox2mm: 1-based [i,j,k] row-vector → world mm row-vector
            vox2mm = @(v) ([v, 1] * Txfm);

            % ---- Sort chanNames (same logic as namingGUI) ----
            sortedChanNames = self.chanNames;
            if ~isempty(sortedChanNames)
                toks = regexp(sortedChanNames, '^(.*)\s+(\d+)$', 'tokens', 'once');
                valid_tok = ~cellfun(@isempty, toks);
                keys = sortedChanNames;
                keys(valid_tok) = cellfun( ...
                    @(t) sprintf('%s_%05d', t{1}, str2double(t{2})), ...
                    toks(valid_tok), 'UniformOutput', false);
                [~, si] = sort(keys);
                sortedChanNames = sortedChanNames(si);
            end

            % ---- Mutable state shared across GUI and restart loops ----
            markers = struct('chanName',{},'type',{},'x',{},'y',{},'z',{});

            % ---- Main loop: slicer → 3-D review → restart or proceed ----
            markersFS = [];   % FreeSurfer RAS coords; set when user proceeds
            while true
                markers = runSlicerGUI(markers, sortedChanNames);
                if isempty(markers) && nargout == 0
                    error('electrodeLocalizer:cancelled', ...
                        '[manualLocalize] Cancelled — no markers placed.');
                end
                % Transform CT NIfTI RAS → scanner RAS, then snap subdural
                % contacts to the pial-outer-smoothed surface before review.
                markersFS = ctRasToFS(markers);
                markersFS = snapToPia(markersFS);
                action = run3DReview(markersFS);
                if strcmp(action, 'proceed')
                    break;
                elseif strcmp(action, 'quit')
                    error('electrodeLocalizer:userQuit', ...
                        '[manualLocalize] Aborted by user in 3-D review.');
                end
                % 'restart' → loop back to slicer; markers stay in CT RAS
            end

            % ---- Build self.leads (FreeSurfer RAS coords) ----
            N = numel(markersFS);
            chanNames_out = {markersFS.chanName}';
            types_out     = {markersFS.type}';
            xyz = [[markersFS.x]', [markersFS.y]', [markersFS.z]'];
            self.clusters = xyz;
            self.leads = table(chanNames_out, xyz(:,1), xyz(:,2), xyz(:,3), types_out, ...
                'VariableNames', {'chanName','x','y','z','type'});

            clustFile = fullfile(self.locDirs.ct_1_xfm, 'clusters_mr.mat');
            clusters_mm = xyz;                                       %#ok<NASGU>
            save(clustFile, 'clusters_mm');
            fprintf('[manualLocalize] %d contacts saved.\n', N);

            % ==============================================================
            % Nested: CT NIfTI RAS → FreeSurfer RAS transform
            % ==============================================================
            function mkrsFS = ctRasToFS(mkrs)
                % Converts marker coords (stored as NIfTI world mm, RAS) to
                % FreeSurfer/scanner RAS (same space as SUMA pial GIFTIs).
                %
                % Pipeline (bypasses xfm_leads/ConvertSurface entirely):
                %   1. CT NIfTI RAS → CT NIfTI 0-idx voxel  (via inv(Txfm))
                %   2. CT NIfTI voxel → CT BRIK voxel        (axis flips for RAI BRIK)
                %   3. CT BRIK voxel → CT AFNI DICOM mm      (scale + origin)
                %   4. CT AFNI DICOM → MR AFNI DICOM         (inv(aff12.1D))
                %   5. MR AFNI DICOM → MR BRIK voxel         (reverse scale)
                %   6. MR BRIK voxel → mr_pre.nii voxel      (all axes flipped: RAI→LPS)
                %   7. mr_pre.nii voxel → NIfTI RAS          (mr_pre.nii affine)
                if isempty(mkrs)
                    mkrsFS = mkrs;
                    return;
                end
                nReal = numel(mkrs);

                % ---- Step 1: CT NIfTI RAS → 0-indexed NIfTI voxel ----
                invT = inv(Txfm);
                vox_0 = zeros(nReal, 3);
                for mi = 1:nReal
                    mm_h = [mkrs(mi).x, mkrs(mi).y, mkrs(mi).z, 1];
                    vox_h = mm_h * invT;
                    vox_0(mi,:) = round(vox_h(1:3)) - 1;
                end

                % ---- Step 2: CT NIfTI voxel → CT BRIK voxel (RAI) ----
                % Flip each axis where NIfTI and BRIK directions are opposite.
                % For LPS NIfTI (Txfm diagonal: neg, neg, pos): all 3 axes flip.
                ct_bvox = vox_0;
                if Txfm(1,1) < 0, ct_bvox(:,1) = (nx-1) - ct_bvox(:,1); end
                if Txfm(2,2) < 0, ct_bvox(:,2) = (ny-1) - ct_bvox(:,2); end
                if Txfm(3,3) > 0, ct_bvox(:,3) = (nz-1) - ct_bvox(:,3); end

                % ---- Step 3: CT BRIK voxel → CT AFNI DICOM ----
                [~, ct_hdr] = unix(sprintf('"%s/3dinfo" -di -dj -dk -o3 "%s+orig" 2>/dev/null', ...
                    self.afniBin, fullfile(workDir,'ct_implant')));
                ct_v = sscanf(ct_hdr, '%f');   % [di dj dk ox oy oz]
                ct_vs  = ct_v(1:3)';
                ct_org = ct_v(4:6)';
                ct_dicom = ct_bvox .* ct_vs + ct_org;

                % ---- Step 4: CT AFNI DICOM → MR AFNI DICOM via inv(aff12.1D) ----
                fid    = fopen(aff1D, 'r');
                M_vals = fscanf(fid, '%f', 12);
                fclose(fid);
                M_raw  = reshape(M_vals(:), 4, 3)';   % 3×4 (column-major fill → transpose)
                M4     = [M_raw; 0 0 0 1];
                invM4  = inv(M4);
                ct_h     = [ct_dicom, ones(nReal,1)]';  % 4 × nReal
                mr_dicom = (invM4 * ct_h)';             % nReal × 4
                mr_dicom = mr_dicom(:,1:3);

                mrBrikDir   = workDir;
                mrNiiForXfm = fullfile(workDir, 'mr_pre.nii');

                % ---- Step 5: MR AFNI DICOM → MR BRIK voxel ----
                [~, mr_hdr] = unix(sprintf('"%s/3dinfo" -ni -nj -nk -di -dj -dk -o3 "%s" 2>/dev/null', ...
                    self.afniBin, fullfile(mrBrikDir,'mr_pre_do+orig')));
                mr_v       = sscanf(mr_hdr, '%f');   % [ni nj nk di dj dk ox oy oz]
                mr_brik_dims = mr_v(1:3)';
                mr_vs        = mr_v(4:6)';
                mr_org       = mr_v(7:9)';
                mr_bvox = round((mr_dicom - mr_org) ./ mr_vs);

                % ---- Step 6: MR BRIK voxel → mr_pre.nii voxel ----
                % RAI BRIK and mr_pre.nii (LPS) have all 3 axes reversed.
                % Use mr_brik_dims (not mr_nii dims) because the de-obliqued
                % BRIK has more z-slices (278) than the original NIfTI (250).
                mr_nii_vox = (mr_brik_dims - 1) - mr_bvox;

                % ---- Step 7: mr_pre.nii voxel → NIfTI RAS ----
                mr_nii_info = niftiinfo(mrNiiForXfm);
                mr_T = mr_nii_info.Transform.T;  % 4×4 MATLAB affine (1-based row-vec → world)
                mr_vox_1 = mr_nii_vox + 1;       % 0-indexed → 1-indexed
                mr_ras_h = [mr_vox_1, ones(nReal,1)] * mr_T;  % nReal × 4
                xyzFS = mr_ras_h(:, 1:3);

                mkrsFS = mkrs;
                for mi = 1:nReal
                    mkrsFS(mi).x = xyzFS(mi,1);
                    mkrsFS(mi).y = xyzFS(mi,2);
                    mkrsFS(mi).z = xyzFS(mi,3);
                end
            end

            % ==============================================================
            % Nested: snap subdural contacts to pial-outer-smoothed surface
            % ==============================================================
            function mkrsOut = snapToPia(mkrsIn)
                mkrsOut = mkrsIn;
                if isempty(mkrsIn), return; end

                sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');
                lhFile  = fullfile(sumaDir, 'lh.pial-outer-smoothed.gii');
                rhFile  = fullfile(sumaDir, 'rh.pial-outer-smoothed.gii');
                if ~exist(lhFile,'file') || ~exist(rhFile,'file')
                    return;   % surfaces not yet available — skip silently
                end
                lhVerts = gifti(lhFile).vertices;
                rhVerts = gifti(rhFile).vertices;

                for mi = 1:numel(mkrsIn)
                    if ~strcmp(mkrsIn(mi).type, 'subdural'), continue; end
                    if mkrsIn(mi).x < 0
                        verts = lhVerts;
                    else
                        verts = rhVerts;
                    end
                    dists = sum((verts - [mkrsIn(mi).x, mkrsIn(mi).y, mkrsIn(mi).z]).^2, 2);
                    [~, idx] = min(dists);
                    mkrsOut(mi).x = verts(idx,1);
                    mkrsOut(mi).y = verts(idx,2);
                    mkrsOut(mi).z = verts(idx,3);
                end
            end

            % ==============================================================
            % Nested: CT slicer GUI
            % ==============================================================
            function markersOut = runSlicerGUI(markersIn, sortedChanNames)

                % ---- State ----
                curVox = [round(nx/2), round(ny/2), round(nz/2)];
                wW = 3000;  wL = 700;
                mode = 'scroll';   % 'scroll' | 'zoom'
                dragStart    = [];  % axes data coords at drag start
                dragStartFig = [];  % figure pixel coords at drag start
                dragLimX     = [];
                dragLimY     = [];
                dragAxes     = [];
                isDragging   = false;
                pixdim       = info.PixelDimensions;   % [dx dy dz] mm/vox
                markersOut   = markersIn;
                pendingPos   = [];   % voxel [ix,iy,iz] awaiting name entry

                % ---- Orientation flip flags (computed once, shared by all callbacks) ----
                % These determine whether the displayed image is flipud relative to
                % voxel index order.  Any callback that maps display coords ↔ voxels
                % must apply the same transforms used in refreshAll.
                colX_ = Txfm(1:3,1);  colY_ = Txfm(1:3,2);  colZ_ = Txfm(1:3,3);
                ax_xFlip  = colX_(1) > 0;   % XDir='reverse' for axial/coronal
                ax_yFlip  = colY_(2) < 0;   % flipud on axial y (dim2)
                cor_xFlip = colX_(1) > 0;
                cor_yTop  = colZ_(3) < 0;   % flipud on coronal/sagittal y (dim3)
                sag_xFlip = colY_(2) < 0;   % XDir='reverse' for sagittal x

                % ---- Color scheme ----
                bg  = [0.12 0.12 0.12];
                bg2 = [0.18 0.18 0.18];
                bg3 = [0.22 0.22 0.22];
                fg  = [0.92 0.92 0.92];
                acc = [0  0.9 0.9];    % cyan accent

                % ---- Figure ----
                fig = figure('Name', sprintf('CT Slicer — %s', self.subj), ...
                    'NumberTitle','off','Color',bg, ...
                    'Position',[30 30 1500 870], ...
                    'KeyPressFcn',@cbKey, ...
                    'WindowScrollWheelFcn',@cbScroll, ...
                    'WindowButtonDownFcn',@cbDown, ...
                    'WindowButtonMotionFcn',@cbMotion, ...
                    'WindowButtonUpFcn',@cbUp, ...
                    'CloseRequestFcn',@cbClose);

                % ---- Toolbar strip (top, left 75%) ----
                tbH  = 0.048;
                tb1Y = 0.950;           % top row: mode buttons + hint
                tb2Y = tb1Y - tbH - 0.006;  % bottom row: W/L controls + Reset View

                % --- Row 1: mode + hint ---
                uicontrol(fig,'Style','text','String','Mode:', ...
                    'Units','normalized','Position',[0.01 tb1Y 0.04 tbH], ...
                    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',10);
                btnNormal = uicontrol(fig,'Style','togglebutton','String','Scroll', ...
                    'Units','normalized','Position',[0.05 tb1Y 0.07 tbH], ...
                    'BackgroundColor',[0.25 0.55 0.25],'ForegroundColor','w', ...
                    'FontSize',10,'Value',1,'Callback',@(~,~) setMode('scroll'));
                btnZoom = uicontrol(fig,'Style','togglebutton','String','Zoom', ...
                    'Units','normalized','Position',[0.13 tb1Y 0.06 tbH], ...
                    'BackgroundColor',bg2,'ForegroundColor',fg, ...
                    'FontSize',10,'Value',0,'Callback',@(~,~) setMode('zoom'));
                uicontrol(fig,'Style','text', ...
                    'String','Drag to pan  |  Click to place crosshair  |  Scroll: slice (scroll mode) / zoom (zoom mode)', ...
                    'Units','normalized','Position',[0.20 tb1Y 0.56 tbH], ...
                    'BackgroundColor',bg,'ForegroundColor',[0.5 0.5 0.5],'FontSize',9, ...
                    'HorizontalAlignment','left');

                % --- Row 2: W/L + Reset View ---
                uicontrol(fig,'Style','text','String','W:', ...
                    'Units','normalized','Position',[0.01 tb2Y 0.025 tbH], ...
                    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',10);
                hW = uicontrol(fig,'Style','edit','String',num2str(wW), ...
                    'Units','normalized','Position',[0.038 tb2Y 0.060 tbH], ...
                    'BackgroundColor',bg3,'ForegroundColor',fg, ...
                    'FontSize',10,'Callback',@cbWL);
                uicontrol(fig,'Style','text','String','L:', ...
                    'Units','normalized','Position',[0.103 tb2Y 0.025 tbH], ...
                    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',10);
                hL = uicontrol(fig,'Style','edit','String',num2str(wL), ...
                    'Units','normalized','Position',[0.131 tb2Y 0.060 tbH], ...
                    'BackgroundColor',bg3,'ForegroundColor',fg, ...
                    'FontSize',10,'Callback',@cbWL);
                uicontrol(fig,'Style','pushbutton','String','Reset View', ...
                    'Units','normalized','Position',[0.20 tb2Y 0.08 tbH], ...
                    'BackgroundColor',bg2,'ForegroundColor',fg, ...
                    'FontSize',10,'Callback',@cbResetView);

                % ---- Slice panels ----
                % Make panels square in pixels so DataAspectRatio fills them
                % without black bars for isotropic CT data.
                figPos = get(fig,'Position');   % [left bottom width height] px
                panW = 0.235;
                panH = panW * figPos(3) / figPos(4);  % same pixel size → square
                panY = (tb2Y - panH) / 2 + 0.005;    % centered between bottom edge and row-2 toolbar
                axAx  = axes('Parent',fig,'Position',[0.005 panY panW panH], ...
                    'Color','k','XColor','none','YColor','none'); hold(axAx,'on');
                axCor = axes('Parent',fig,'Position',[0.248 panY panW panH], ...
                    'Color','k','XColor','none','YColor','none'); hold(axCor,'on');
                axSag = axes('Parent',fig,'Position',[0.491 panY panW panH], ...
                    'Color','k','XColor','none','YColor','none'); hold(axSag,'on');

                % Image handles (replaced on each refresh)
                hImAx  = []; hImCor = []; hImSag = [];

                % Crosshair handles
                hXHax  = [];  hXHcor = [];  hXHsag = [];

                % Marker scatter handles
                hMkAx  = [];  hMkCor = [];  hMkSag = [];

                % Slice index labels
                hTitleAx  = title(axAx,  '','Color',fg,'FontSize',9);
                hTitleCor = title(axCor, '','Color',fg,'FontSize',9);
                hTitleSag = title(axSag, '','Color',fg,'FontSize',9);

                % Status bar
                hStatus = uicontrol(fig,'Style','text', ...
                    'Units','normalized','Position',[0.005 0.003 0.73 0.038], ...
                    'BackgroundColor',bg,'ForegroundColor',[0.5 0.8 1.0], ...
                    'FontSize',10,'HorizontalAlignment','left','String','Click on a contact location');

                % ---- Right panel ----
                rx = 0.742;  rw = 0.253;
                uicontrol(fig,'Style','text','String','CT Slicer', ...
                    'Units','normalized','Position',[rx 0.94 rw 0.05], ...
                    'BackgroundColor',bg,'ForegroundColor',acc, ...
                    'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');

                uicontrol(fig,'Style','text','String','Channel name:', ...
                    'Units','normalized','Position',[rx 0.89 rw 0.04], ...
                    'BackgroundColor',bg,'ForegroundColor',fg, ...
                    'FontSize',10,'HorizontalAlignment','left');
                hList = uicontrol(fig,'Style','listbox', ...
                    'Units','normalized','Position',[rx 0.68 rw 0.21], ...
                    'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',11,'String',{});
                uicontrol(fig,'Style','text','String','Or type name:', ...
                    'Units','normalized','Position',[rx 0.635 rw 0.04], ...
                    'BackgroundColor',bg,'ForegroundColor',fg, ...
                    'FontSize',10,'HorizontalAlignment','left');
                hEdit = uicontrol(fig,'Style','edit', ...
                    'Units','normalized','Position',[rx 0.585 rw 0.048], ...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',fg, ...
                    'FontSize',11,'HorizontalAlignment','left','String','');
                uicontrol(fig,'Style','text','String','Electrode type:', ...
                    'Units','normalized','Position',[rx 0.54 rw 0.04], ...
                    'BackgroundColor',bg,'ForegroundColor',fg, ...
                    'FontSize',10,'HorizontalAlignment','left');
                hType = uicontrol(fig,'Style','popupmenu', ...
                    'Units','normalized','Position',[rx 0.49 rw 0.048], ...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',fg, ...
                    'FontSize',11,'String',{'depth','subdural'},'Value',1);
                hPlaceBtn = uicontrol(fig,'Style','pushbutton', ...
                    'String','Place Marker', ...
                    'Units','normalized','Position',[rx 0.42 rw 0.063], ...
                    'BackgroundColor',[0.18 0.48 0.18],'ForegroundColor','w', ...
                    'FontSize',12,'FontWeight','bold','Enable','off','Callback',@cbPlace);

                uicontrol(fig,'Style','text','String','Placed markers:', ...
                    'Units','normalized','Position',[rx 0.37 rw 0.04], ...
                    'BackgroundColor',bg,'ForegroundColor',fg, ...
                    'FontSize',10,'HorizontalAlignment','left');
                hMkList = uicontrol(fig,'Style','listbox', ...
                    'Units','normalized','Position',[rx 0.19 rw 0.18], ...
                    'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',10,'String',{});
                uicontrol(fig,'Style','pushbutton','String','Remove Selected', ...
                    'Units','normalized','Position',[rx 0.13 rw 0.056], ...
                    'BackgroundColor',[0.48 0.18 0.18],'ForegroundColor','w', ...
                    'FontSize',10,'Callback',@cbRemove);
                uicontrol(fig,'Style','pushbutton','String','Quit', ...
                    'Units','normalized','Position',[rx 0.02 rw*0.44 0.09], ...
                    'BackgroundColor',[0.35 0.15 0.15],'ForegroundColor','w', ...
                    'FontSize',11,'Callback',@cbQuit);
                uicontrol(fig,'Style','pushbutton','String','Done', ...
                    'Units','normalized','Position',[rx+rw*0.56 0.02 rw*0.44 0.09], ...
                    'BackgroundColor',[0.20 0.38 0.58],'ForegroundColor','w', ...
                    'FontSize',11,'FontWeight','bold','Callback',@cbDone);

                refreshAll();
                waitfor(fig);

                % ---- Helpers ----
                function setMode(m)
                    mode = m;
                    set(btnNormal,'Value', strcmp(m,'scroll'), ...
                        'BackgroundColor', ternary(strcmp(m,'scroll'),[0.25 0.55 0.25],bg2));
                    set(btnZoom,  'Value', strcmp(m,'zoom'), ...
                        'BackgroundColor', ternary(strcmp(m,'zoom'), [0.10 0.35 0.55],bg2));
                end

                function v = ternary(cond, a, b)
                    if cond, v = a; else, v = b; end
                end

                function addOrientLabels(ax, leftLbl, rightLbl, topLbl, botLbl)
                    % Place orientation labels at the visual edges of ax.
                    % Uses DATA units to avoid any ambiguity with reversed axes:
                    % XDir='reverse' means xl(2) is the visual-left edge.
                    prevLbls = findobj(ax,'Tag','orientLabel');
                    if ~isempty(prevLbls), delete(prevLbls); end
                    xl = get(ax,'XLim');  yl = get(ax,'YLim');
                    xrev = strcmp(get(ax,'XDir'),'reverse');
                    yrev = strcmp(get(ax,'YDir'),'reverse');
                    % Visual-edge data coordinates.
                    xL = xl(1 + xrev);   % xl(2) if reversed (high x = visual left)
                    xR = xl(2 - xrev);   % xl(1) if reversed (low  x = visual right)
                    yT = yl(2 - yrev);   % yl(2) if normal   (high y = visual top)
                    yB = yl(1 + yrev);   % yl(1) if normal   (low  y = visual bot)
                    % Small inset (move away from the edge, toward center).
                    dx = 0.03*(xl(2)-xl(1)) * (1 - 2*xrev);  % +dx inward if normal, -dx if rev
                    dy = 0.03*(yl(2)-yl(1)) * (1 - 2*yrev);
                    txtArgs = {'Units','data','FontSize',10,'FontWeight','bold', ...
                               'Color',[1 1 0],'HitTest','off','Tag','orientLabel', ...
                               'Clipping','off'};
                    text(ax, xL+dx, mean(yl), leftLbl,  txtArgs{:}, ...
                        'HorizontalAlignment','left',  'VerticalAlignment','middle');
                    text(ax, xR-dx, mean(yl), rightLbl, txtArgs{:}, ...
                        'HorizontalAlignment','right', 'VerticalAlignment','middle');
                    text(ax, mean(xl), yT-dy, topLbl,   txtArgs{:}, ...
                        'HorizontalAlignment','center','VerticalAlignment','top');
                    text(ax, mean(xl), yB+dy, botLbl,   txtArgs{:}, ...
                        'HorizontalAlignment','center','VerticalAlignment','bottom');
                end

                function ax = hitAxes(pt)
                    % pt = figure CurrentPoint [x y] pixels
                    figPos = get(fig,'Position');
                    np = pt ./ figPos(3:4);
                    ax = [];
                    for aa = [axAx, axCor, axSag]
                        apos = get(aa,'Position');
                        if np(1) >= apos(1) && np(1) <= apos(1)+apos(3) && ...
                           np(2) >= apos(2) && np(2) <= apos(2)+apos(4)
                            ax = aa;
                            return;
                        end
                    end
                end

                function [dp, valid] = axesDataPoint(ax)
                    cp = get(ax,'CurrentPoint');
                    dp = cp(1,1:2);
                    xl = get(ax,'XLim');  yl = get(ax,'YLim');
                    valid = dp(1)>=xl(1) && dp(1)<=xl(2) && dp(2)>=yl(1) && dp(2)<=yl(2);
                end

                function refreshAll()
                    if ~ishandle(fig), return; end
                    clim_lo = wL - wW/2;
                    clim_hi = wL + wW/2;

                    % Flip flags are computed once in the outer scope (ax_xFlip, etc.)

                    % ---- Axial ----
                    sl_ax = vol(:,:,curVox(3))';
                    if ax_yFlip, sl_ax = flipud(sl_ax); end
                    if isempty(hImAx) || ~ishandle(hImAx)
                        hImAx = imagesc(axAx, sl_ax, [clim_lo clim_hi]);
                        colormap(axAx, gray);
                        axis(axAx,'normal');
                        set(axAx,'XLim',[0.5 nx+0.5],'YLim',[0.5 ny+0.5],'YDir','normal');
                        if ax_xFlip, set(axAx,'XDir','reverse'); end
                        addOrientLabels(axAx, 'R','L','A','P');
                    else
                        set(hImAx,'CData',sl_ax);
                        clim(axAx,[clim_lo clim_hi]);
                    end
                    set(hTitleAx,'String',sprintf('AXIAL   z=%d/%d',curVox(3),nz));

                    % ---- Coronal ----
                    sl_cor = squeeze(vol(:,curVox(2),:))';
                    if cor_yTop, sl_cor = flipud(sl_cor); end
                    if isempty(hImCor) || ~ishandle(hImCor)
                        hImCor = imagesc(axCor, sl_cor, [clim_lo clim_hi]);
                        colormap(axCor, gray);
                        axis(axCor,'normal');
                        set(axCor,'XLim',[0.5 nx+0.5],'YLim',[0.5 nz+0.5],'YDir','normal');
                        if cor_xFlip, set(axCor,'XDir','reverse'); end
                        addOrientLabels(axCor, 'R','L','S','I');
                    else
                        set(hImCor,'CData',sl_cor);
                        clim(axCor,[clim_lo clim_hi]);
                    end
                    set(hTitleCor,'String',sprintf('CORONAL   y=%d/%d',curVox(2),ny));

                    % ---- Sagittal ----
                    sl_sag = squeeze(vol(curVox(1),:,:))';
                    if cor_yTop, sl_sag = flipud(sl_sag); end   % same z-flip as coronal
                    if isempty(hImSag) || ~ishandle(hImSag)
                        hImSag = imagesc(axSag, sl_sag, [clim_lo clim_hi]);
                        colormap(axSag, gray);
                        axis(axSag,'normal');
                        set(axSag,'XLim',[0.5 ny+0.5],'YLim',[0.5 nz+0.5],'YDir','normal');
                        if sag_xFlip, set(axSag,'XDir','reverse'); end
                        addOrientLabels(axSag, 'P','A','S','I');
                    else
                        set(hImSag,'CData',sl_sag);
                        clim(axSag,[clim_lo clim_hi]);
                    end
                    set(hTitleSag,'String',sprintf('SAGITTAL   x=%d/%d',curVox(1),nx));

                    refreshCrosshairs();
                    refreshMarkers();

                    mm = vox2mm(curVox);
                    set(hStatus,'String', ...
                        sprintf('  vox [%d, %d, %d]    x=%.1f  y=%.1f  z=%.1f mm', ...
                        curVox(1),curVox(2),curVox(3), mm(1),mm(2),mm(3)));

                    % Refresh available channel names list
                    if ~isempty(sortedChanNames)
                        used = {markersOut.chanName};
                        rem  = sortedChanNames(~ismember(sortedChanNames, used));
                        set(hList,'String',rem(:), ...
                            'Value',max(1,min(get(hList,'Value'),numel(rem))));
                    end
                end

                function refreshCrosshairs()
                    if ~ishandle(fig), return; end
                    % Axial crosshair: vertical at ix, horizontal at iy.
                    % When ax_yFlip=true, the image was flipud before display,
                    % so voxel iy appears at display row (ny - iy + 1).
                    iy_disp = ternary(ax_yFlip, ny - curVox(2) + 1, curVox(2));
                    if ~isempty(hXHax) && all(ishandle(hXHax))
                        delete(hXHax);
                    end
                    hXHax(1) = plot(axAx, [curVox(1) curVox(1)], [0.5 ny+0.5], ...
                        'Color',[1 0.3 0.3],'LineWidth',0.8,'HitTest','off');
                    hXHax(2) = plot(axAx, [0.5 nx+0.5], [iy_disp iy_disp], ...
                        'Color',[1 0.3 0.3],'LineWidth',0.8,'HitTest','off');

                    % Coronal crosshair: vertical at ix, horizontal at (flipped) iz
                    % When cor_yTop=true the image was flipud, so voxel iz maps to display row (nz - iz + 1).
                    iz_cor = ternary(cor_yTop, nz - curVox(3) + 1, curVox(3));
                    if ~isempty(hXHcor) && all(ishandle(hXHcor))
                        delete(hXHcor);
                    end
                    hXHcor(1) = plot(axCor, [curVox(1) curVox(1)], [0.5 nz+0.5], ...
                        'Color',[1 0.3 0.3],'LineWidth',0.8,'HitTest','off');
                    hXHcor(2) = plot(axCor, [0.5 nx+0.5], [iz_cor iz_cor], ...
                        'Color',[1 0.3 0.3],'LineWidth',0.8,'HitTest','off');

                    % Sagittal crosshair: vertical at iy, horizontal at (flipped) iz
                    if ~isempty(hXHsag) && all(ishandle(hXHsag))
                        delete(hXHsag);
                    end
                    hXHsag(1) = plot(axSag, [curVox(2) curVox(2)], [0.5 nz+0.5], ...
                        'Color',[1 0.3 0.3],'LineWidth',0.8,'HitTest','off');
                    hXHsag(2) = plot(axSag, [0.5 ny+0.5], [iz_cor iz_cor], ...
                        'Color',[1 0.3 0.3],'LineWidth',0.8,'HitTest','off');
                end

                function refreshMarkers()
                    if ~ishandle(fig), return; end
                    if ~isempty(hMkAx)  && all(ishandle(hMkAx)),  delete(hMkAx);  end
                    if ~isempty(hMkCor) && all(ishandle(hMkCor)), delete(hMkCor); end
                    if ~isempty(hMkSag) && all(ishandle(hMkSag)), delete(hMkSag); end
                    hMkAx = []; hMkCor = []; hMkSag = [];

                    if isempty(markersOut), return; end

                    % Convert marker world mm back to voxel for display.
                    % vox2mm: vox_h * Txfm = mm_h  →  inverse: mm_h * inv(Txfm) = vox_h
                    invT = inv(Txfm);
                    for mi = 1:numel(markersOut)
                        mm_h = [markersOut(mi).x, markersOut(mi).y, markersOut(mi).z, 1];
                        vx   = round(mm_h * invT);   % no transpose — row-vector convention
                        vx   = vx(1:3);

                        mkCol  = [0.15 0.35 0.85];   % dark blue dot
                        txtCol = [0.55 0.75 1.00];   % lighter blue text (readable on dark CT)

                        % Axial: only show if this marker is on the current axial slice (z).
                        if abs(vx(3) - curVox(3)) <= 1
                            vy_disp = ternary(ax_yFlip, ny - vx(2) + 1, vx(2));
                            hMkAx(end+1) = scatter(axAx, vx(1), vy_disp, 80, mkCol, ...
                                'filled','HitTest','off','MarkerEdgeColor','w','LineWidth',0.5); %#ok<AGROW>
                            hMkAx(end+1) = text(axAx, vx(1)+2, vy_disp, ...
                                markersOut(mi).chanName,'Color',txtCol,'FontSize',8, ...
                                'FontWeight','bold','HitTest','off'); %#ok<AGROW>
                        end

                        % Coronal: only show if on the current coronal slice (y).
                        iz_c = ternary(cor_yTop, nz - vx(3) + 1, vx(3));
                        if abs(vx(2) - curVox(2)) <= 1
                            hMkCor(end+1) = scatter(axCor, vx(1), iz_c, 80, mkCol, ...
                                'filled','HitTest','off','MarkerEdgeColor','w','LineWidth',0.5); %#ok<AGROW>
                            hMkCor(end+1) = text(axCor, vx(1)+2, iz_c, ...
                                markersOut(mi).chanName,'Color',txtCol,'FontSize',8, ...
                                'FontWeight','bold','HitTest','off'); %#ok<AGROW>
                        end

                        % Sagittal: only show if on the current sagittal slice (x).
                        if abs(vx(1) - curVox(1)) <= 1
                            hMkSag(end+1) = scatter(axSag, vx(2), iz_c, 80, mkCol, ...
                                'filled','HitTest','off','MarkerEdgeColor','w','LineWidth',0.5); %#ok<AGROW>
                            hMkSag(end+1) = text(axSag, vx(2)+2, iz_c, ...
                                markersOut(mi).chanName,'Color',txtCol,'FontSize',8, ...
                                'FontWeight','bold','HitTest','off'); %#ok<AGROW>
                        end
                    end

                    % Update marker list on right panel
                    strs = arrayfun(@(m) sprintf('%s  [%s]', m.chanName, m.type), ...
                        markersOut, 'UniformOutput', false);
                    set(hMkList,'String',strs,'Value',min(get(hMkList,'Value'),numel(strs)));
                end

                function setPendingPos(vox)
                    pendingPos = vox;
                    curVox = vox;
                    set(hPlaceBtn,'Enable','on');
                    recenterAllOnCurVox();
                    refreshAll();
                end

                function recenterAllOnCurVox()
                    % Pan each panel to keep curVox centered, preserving zoom level.
                    vy_disp = ternary(ax_yFlip, ny - curVox(2) + 1, curVox(2));
                    iz_disp = ternary(cor_yTop, nz - curVox(3) + 1, curVox(3));
                    centerAxis(axAx,  curVox(1), vy_disp);
                    centerAxis(axCor, curVox(1), iz_disp);
                    centerAxis(axSag, curVox(2), iz_disp);
                    addOrientLabels(axAx,  'R','L','A','P');
                    addOrientLabels(axCor, 'R','L','S','I');
                    addOrientLabels(axSag, 'P','A','S','I');
                end

                function centerAxis(ax, cx, cy)
                    xl = get(ax,'XLim');  yl = get(ax,'YLim');
                    hw = (xl(2)-xl(1))/2;  hh = (yl(2)-yl(1))/2;
                    set(ax,'XLim',[cx-hw, cx+hw], 'YLim',[cy-hh, cy+hh]);
                end

                % ---- Window-level callbacks ----
                function cbWL(~,~)
                    wW_new = str2double(get(hW,'String'));
                    wL_new = str2double(get(hL,'String'));
                    if ~isnan(wW_new) && wW_new > 0, wW = wW_new; end
                    if ~isnan(wL_new), wL = wL_new; end
                    refreshAll();
                end

                % ---- Scroll ----
                function cbScroll(~, evt)
                    ax = hitAxes(get(fig,'CurrentPoint'));
                    delta = -evt.VerticalScrollCount;
                    if isempty(ax)
                        return
                    elseif strcmp(mode,'zoom') && ~isempty(ax)
                        factor = 1 + 0.15*delta;
                        % Use cursor position as center for the hovered panel;
                        % use curVox as center for the other two panels.
                        cp = get(ax,'CurrentPoint');
                        ax_cx = cp(1,1);  ax_cy = cp(1,2);
                        vy_disp = ternary(ax_yFlip, ny - curVox(2) + 1, curVox(2));
                        iz_disp = ternary(cor_yTop, nz - curVox(3) + 1, curVox(3));
                        % Axial center: (vox_x, display_y)
                        if isequal(ax, axAx)
                            zoomAxes(axAx, factor, ax_cx, ax_cy);
                        else
                            zoomAxes(axAx, factor, curVox(1), vy_disp);
                        end
                        % Coronal center: (vox_x, display_z)
                        if isequal(ax, axCor)
                            zoomAxes(axCor, factor, ax_cx, ax_cy);
                        else
                            zoomAxes(axCor, factor, curVox(1), iz_disp);
                        end
                        % Sagittal center: (vox_y, display_z)
                        if isequal(ax, axSag)
                            zoomAxes(axSag, factor, ax_cx, ax_cy);
                        else
                            zoomAxes(axSag, factor, curVox(2), iz_disp);
                        end
                        return;
                    end
                    if isequal(ax, axAx)
                        curVox(3) = max(1, min(nz, curVox(3) + delta));
                    elseif isequal(ax, axCor)
                        curVox(2) = max(1, min(ny, curVox(2) + delta));
                    elseif isequal(ax, axSag)
                        curVox(1) = max(1, min(nx, curVox(1) + delta));
                    end
                    refreshAll();
                end

                function zoomAxes(ax, factor, cx, cy)
                    xl  = get(ax,'XLim');  yl = get(ax,'YLim');
                    nxl = cx + (xl - cx) / factor;
                    nyl = cy + (yl - cy) / factor;
                    set(ax,'XLim',nxl,'YLim',nyl);
                    % Reposition orientation labels to current view edges
                    if isequal(ax, axAx)
                        addOrientLabels(axAx,  'R','L','A','P');
                    elseif isequal(ax, axCor)
                        addOrientLabels(axCor, 'R','L','S','I');
                    elseif isequal(ax, axSag)
                        addOrientLabels(axSag, 'P','A','S','I');
                    end
                end

                % ---- Drag to pan / click to crosshair ----
                function cbDown(~,~)
                    ax = hitAxes(get(fig,'CurrentPoint'));
                    if isempty(ax), return; end
                    dragAxes     = ax;
                    cp           = get(ax,'CurrentPoint');
                    dragStart    = cp(1,1:2);
                    dragLimX     = get(ax,'XLim');
                    dragLimY     = get(ax,'YLim');
                    dragStartFig = get(fig,'CurrentPoint');
                    isDragging   = false;
                end
                function cbMotion(~,~)
                    if isempty(dragStart) || isempty(dragAxes), return; end
                    figPt = get(fig,'CurrentPoint');
                    if norm(figPt - dragStartFig) > 4
                        isDragging = true;
                        cp    = get(dragAxes,'CurrentPoint');
                        delta = cp(1,1:2) - dragStart;
                        set(dragAxes,'XLim', dragLimX - delta(1), ...
                                     'YLim', dragLimY - delta(2));
                        % Reposition orientation labels to current view edges
                        if isequal(dragAxes, axAx)
                            addOrientLabels(axAx,  'R','L','A','P');
                        elseif isequal(dragAxes, axCor)
                            addOrientLabels(axCor, 'R','L','S','I');
                        elseif isequal(dragAxes, axSag)
                            addOrientLabels(axSag, 'P','A','S','I');
                        end
                    end
                end
                function cbUp(~,~)
                    if ~isDragging && ~isempty(dragAxes)
                        [dp, valid] = axesDataPoint(dragAxes);
                        if valid
                            if isequal(dragAxes, axAx)
                                ix = max(1,min(nx, round(dp(1))));
                                % dp(2) is the display row (after flipud if ax_yFlip).
                                % Invert to get original voxel index.
                                iy_d = max(1,min(ny, round(dp(2))));
                                iy   = ternary(ax_yFlip, ny - iy_d + 1, iy_d);
                                setPendingPos([ix, iy, curVox(3)]);
                            elseif isequal(dragAxes, axCor)
                                ix   = max(1,min(nx, round(dp(1))));
                                iz_f = max(1,min(nz, round(dp(2))));
                                % iz_f is display row; invert only if flipud was applied
                                iz   = ternary(cor_yTop, nz - iz_f + 1, iz_f);
                                setPendingPos([ix, curVox(2), iz]);
                            elseif isequal(dragAxes, axSag)
                                iy   = max(1,min(ny, round(dp(1))));
                                iz_f = max(1,min(nz, round(dp(2))));
                                iz   = ternary(cor_yTop, nz - iz_f + 1, iz_f);
                                setPendingPos([curVox(1), iy, iz]);
                            end
                        end
                    end
                    dragStart = []; dragAxes = []; dragLimX = []; dragLimY = [];
                    dragStartFig = []; isDragging = false;
                end

                % ---- Reset view ----
                function cbResetView(~,~)
                    % Flip flags are in outer scope (ax_xFlip, cor_xFlip, sag_xFlip, cor_yTop)
                    axis(axAx,'normal');
                    set(axAx,'XLim',[0.5 nx+0.5],'YLim',[0.5 ny+0.5],'YDir','normal');
                    if ax_xFlip,  set(axAx,'XDir','reverse');  else, set(axAx,'XDir','normal');  end
                    axis(axCor,'normal');
                    set(axCor,'XLim',[0.5 nx+0.5],'YLim',[0.5 nz+0.5],'YDir','normal');
                    if cor_xFlip, set(axCor,'XDir','reverse'); else, set(axCor,'XDir','normal'); end
                    axis(axSag,'normal');
                    set(axSag,'XLim',[0.5 ny+0.5],'YLim',[0.5 nz+0.5],'YDir','normal');
                    if sag_xFlip, set(axSag,'XDir','reverse'); else, set(axSag,'XDir','normal'); end
                    addOrientLabels(axAx,  'R','L','A','P');
                    addOrientLabels(axCor, 'R','L','S','I');
                    addOrientLabels(axSag, 'P','A','S','I');
                end

                % ---- Keyboard shortcuts ----
                function cbKey(~, evt)
                    switch evt.Key
                        case 'uparrow'
                            ax = hitAxes(get(fig,'CurrentPoint'));
                            if isequal(ax,axAx)
                                curVox(3) = min(nz, curVox(3)+1);
                            elseif isequal(ax,axCor)
                                curVox(2) = min(ny, curVox(2)+1);
                            elseif isequal(ax,axSag)
                                curVox(1) = min(nx, curVox(1)+1);
                            end
                            refreshAll();
                        case 'downarrow'
                            ax = hitAxes(get(fig,'CurrentPoint'));
                            if isequal(ax,axAx)
                                curVox(3) = max(1, curVox(3)-1);
                            elseif isequal(ax,axCor)
                                curVox(2) = max(1, curVox(2)-1);
                            elseif isequal(ax,axSag)
                                curVox(1) = max(1, curVox(1)-1);
                            end
                            refreshAll();
                    end
                end

                % ---- Place marker ----
                function cbPlace(~,~)
                    if isempty(pendingPos), return; end
                    name = strtrim(get(hEdit,'String'));
                    if isempty(name)
                        listStr = get(hList,'String');
                        if ~isempty(listStr)
                            name = listStr{max(1, get(hList,'Value'))};
                        end
                    end
                    if isempty(name)
                        msgbox('Enter a channel name or select from the list.','','warn');
                        return;
                    end
                    typeStrs = get(hType,'String');
                    typ      = typeStrs{get(hType,'Value')};
                    mm       = vox2mm(pendingPos);
                    markersOut(end+1) = struct('chanName',name,'type',typ, ...
                        'x',mm(1),'y',mm(2),'z',mm(3));
                    set(hEdit,'String','');
                    pendingPos = [];
                    set(hPlaceBtn,'Enable','off');
                    refreshAll();
                end

                % ---- Remove marker ----
                function cbRemove(~,~)
                    if isempty(markersOut), return; end
                    idx = get(hMkList,'Value');
                    if idx < 1 || idx > numel(markersOut), return; end
                    markersOut(idx) = [];
                    set(hMkList,'Value', max(1, idx-1));
                    refreshMarkers();
                end

                % ---- Done / close ----
                function cbDone(~,~)
                    if numel(markersOut) == 0
                        choice = questdlg('No markers placed. Exit anyway?', ...
                            'No Markers','Exit','Cancel','Cancel');
                        if ~strcmp(choice,'Exit'), return; end
                    end
                    delete(fig);
                end
                function cbQuit(~,~)
                    choice = questdlg('Quit localization? All placed markers will be discarded.', ...
                        'Quit','Quit','Cancel','Cancel');
                    if ~strcmp(choice,'Quit'), return; end
                    markersOut = struct('chanName',{},'type',{},'x',{},'y',{},'z',{});
                    delete(fig);
                end
                function cbClose(~,~)
                    delete(fig);
                end

            end  % runSlicerGUI

            % ==============================================================
            % Nested: 3-D review
            % ==============================================================
            function action = run3DReview(mkrs)
                action = 'restart';   % default if window is closed

                if isempty(mkrs)
                    action = 'restart';
                    return;
                end

                sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');
                lhFile  = fullfile(sumaDir, 'lh.pial.gii');
                rhFile  = fullfile(sumaDir, 'rh.pial.gii');

                bg  = [0.12 0.12 0.12];
                fg  = [0.92 0.92 0.92];

                fig3 = figure('Name', sprintf('3-D Review — %s', self.subj), ...
                    'NumberTitle','off','Color',bg, ...
                    'Position',[100 80 1200 820], ...
                    'CloseRequestFcn',@cbQuit);

                ax3 = axes('Parent',fig3,'Position',[0.01 0.01 0.73 0.97], ...
                    'Color','k','XColor','none','YColor','none','ZColor','none');
                hold(ax3,'on'); axis(ax3,'equal','off');
                view(ax3,3); camlight(ax3,'headlight'); material(ax3,'dull');

                if exist(lhFile,'file')==2
                    lhS = gifti(lhFile);
                    patch(ax3,'Faces',lhS.faces,'Vertices',lhS.vertices, ...
                        'FaceColor',[0.75 0.70 0.65],'EdgeColor','none', ...
                        'FaceAlpha',0.35,'PickableParts','none','HitTest','off');
                end
                if exist(rhFile,'file')==2
                    rhS = gifti(rhFile);
                    patch(ax3,'Faces',rhS.faces,'Vertices',rhS.vertices, ...
                        'FaceColor',[0.75 0.70 0.65],'EdgeColor','none', ...
                        'FaceAlpha',0.35,'PickableParts','none','HitTest','off');
                end

                mkCol  = [0.15 0.35 0.85];
                txtCol = [0.55 0.75 1.00];
                for mi = 1:numel(mkrs)
                    scatter3(ax3, mkrs(mi).x, mkrs(mi).y, mkrs(mi).z, ...
                        80, mkCol, 'filled','HitTest','off', ...
                        'MarkerEdgeColor','w','LineWidth',0.5);
                    text(ax3, mkrs(mi).x+1, mkrs(mi).y, mkrs(mi).z, ...
                        mkrs(mi).chanName,'Color',txtCol,'FontSize',8, ...
                        'FontWeight','bold','HitTest','off');
                end
                rotate3d(ax3,'on');

                function v = ternary3(cond,a,b)
                    if cond, v=a; else, v=b; end
                end

                % Right panel
                rx = 0.755;  rw = 0.235;
                uicontrol(fig3,'Style','text','String','Review Placements', ...
                    'Units','normalized','Position',[rx 0.88 rw 0.06], ...
                    'BackgroundColor',bg,'ForegroundColor',[0 0.9 0.9], ...
                    'FontSize',13,'FontWeight','bold','HorizontalAlignment','center');
                uicontrol(fig3,'Style','text', ...
                    'String',sprintf('%d contacts placed.\n\nDrag to rotate.\nVerify positions.\nThen Proceed.', numel(mkrs)), ...
                    'Units','normalized','Position',[rx 0.65 rw 0.22], ...
                    'BackgroundColor',bg,'ForegroundColor',fg, ...
                    'FontSize',11,'HorizontalAlignment','left');
                uicontrol(fig3,'Style','pushbutton','String','Proceed', ...
                    'Units','normalized','Position',[rx 0.45 rw 0.10], ...
                    'BackgroundColor',[0.18 0.48 0.18],'ForegroundColor','w', ...
                    'FontSize',13,'FontWeight','bold','Callback',@cbProceed);
                uicontrol(fig3,'Style','pushbutton','String','Restart Slicer', ...
                    'Units','normalized','Position',[rx 0.32 rw 0.10], ...
                    'BackgroundColor',[0.40 0.20 0.10],'ForegroundColor','w', ...
                    'FontSize',12,'Callback',@cbRestart);
                uicontrol(fig3,'Style','pushbutton','String','Quit', ...
                    'Units','normalized','Position',[rx 0.19 rw 0.10], ...
                    'BackgroundColor',[0.55 0.05 0.05],'ForegroundColor','w', ...
                    'FontSize',12,'FontWeight','bold','Callback',@cbQuit);

                waitfor(fig3);

                function cbProceed(~,~)
                    action = 'proceed';
                    delete(fig3);
                end
                function cbRestart(~,~)
                    action = 'restart';
                    delete(fig3);
                end
                function cbQuit(~,~)
                    action = 'quit';
                    delete(fig3);
                end
            end  % run3DReview

        end  % manualLocalize

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
                % Import only the currently selected file; update its row
                % in-place and stay on the dialog.
                k = get(hList, 'Value');
                if present(k), return; end   % already imported
                if strcmp(names{k}, 'leads.csv')
                    filt = {'*.csv', 'CSV file (*.csv)'};
                else
                    filt = {'*.gii', 'GIFTI surface (*.gii)'};
                end
                [f, d] = uigetfile(filt, sprintf('Select %s', names{k}));
                figure(fig);   % restore focus after uigetfile
                if isequal(f, 0), return; end   % cancelled — leave dialog open
                src = fullfile(d, f);
                try
                    if strcmp(names{k}, 'leads.csv')
                        electrodeLocalizer.validateLeadsCSV(src, self.chanNames);
                    else
                        electrodeLocalizer.validateGifti(src);
                    end
                catch e
                    warndlg(e.message, sprintf('Validation failed — %s', names{k}));
                    return;
                end
                if ~exist(talDir,  'dir'), mkdir(talDir);  end
                if ~exist(sumaDir, 'dir'), mkdir(sumaDir); end
                copyfile(src, dests{k});
                fprintf('[import] %s → %s\n', src, dests{k});
                present(k) = true;
                listStrs{k} = ['[OK]  ' names{k}];
                set(hList, 'String', listStrs);
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
