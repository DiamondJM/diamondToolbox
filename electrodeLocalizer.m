
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
        CT_HU_THRESHOLD   = 1500;  % HU threshold for electrode detection (Hounsfield units)
        AFNI_CLUSTER_RMM  = 1;     % 3dclust: max neighbor distance in mm (rmm)
        AFNI_CLUSTER_VMUL = 10;    % 3dclust: minimum cluster volume in µl (vmul)
        AFNI_BIN_DEFAULT  = fullfile(char(java.lang.System.getProperty('user.home')), 'abin');
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
            % Returns true if tal/leads.csv and SUMA gifti surfaces exist.

            leadsFile  = fullfile(self.rootFolder, self.subj, 'tal', 'leads.csv');
            sumaDir    = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhGii      = fullfile(sumaDir, 'lh.pial-outer-smoothed.gii');
            rhGii      = fullfile(sumaDir, 'rh.pial-outer-smoothed.gii');

            tf = exist(leadsFile, 'file') == 2 && ...
                 exist(lhGii,    'file') == 2 && ...
                 exist(rhGii,    'file') == 2;
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

            % If output files are missing, offer the user a choice before
            % committing to the full pipeline.
            if ~self.isComplete()
                choice = questdlg( ...
                    sprintf('Localization files (leads.csv, lh.pial.gii, rh.pial.gii) not found for %s. How would you like to proceed?', self.subj), ...
                    'electrodeLocalizer', ...
                    'Create these files', 'Import these files', 'Cancel', ...
                    'Create these files');

                if isempty(choice) || strcmp(choice, 'Cancel')
                    fprintf('[electrodeLocalizer] Cancelled.\n');
                    return;
                end

                if strcmp(choice, 'Import these files')
                    self.importExistingLocalization();
                    if self.isComplete()
                        return;   % import covered everything needed
                    end
                    % If still incomplete after import (e.g. user only
                    % imported surfaces but not leads), fall through to
                    % run the remaining pipeline stages.
                end
            end

            self.checkPrerequisites('errorIfMissing', true);
            self.getInputFiles();
            self.runSurface();           % never force — recon-all takes hours
            self.runSuma();              % never force — SUMA takes minutes
            self.coregisterCT('forceNew', forceNew);
            self.detectElectrodes('forceNew', forceNew);
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
                fprintf('Select pre-operative MRI (T1 MPRAGE) — .nii, .nii.gz, or .mgz...\n');
                [f, d] = uigetfile(filter, 'Select pre-op MRI');
                if isequal(f, 0)
                    error('[electrodeLocalizer] MRI is required.');
                end
                self.convertToNii(fullfile(d, f), mrDest);
                fprintf('MRI ready at %s\n', mrDest);
            else
                fprintf('MRI already present: %s\n', mrDest);
            end

            % CT
            if exist(ctDest, 'file') ~= 2
                fprintf('Select post-operative CT — .nii, .nii.gz, or .mgz...\n');
                [f, d] = uigetfile(filter, 'Select post-op CT');
                if isequal(f, 0)
                    error('[electrodeLocalizer] CT is required.');
                end
                self.convertToNii(fullfile(d, f), ctDest);
                fprintf('CT ready at %s\n', ctDest);
            else
                fprintf('CT already present: %s\n', ctDest);
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
            done    = exist(lhPial, 'file') == 2 && exist(rhPial, 'file') == 2;

            if done && ~forceNew
                fprintf('[Stage 3] Surface already exists; skipping recon-all.\n');
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

            % align.sh expects filenames without extension and a work dir.
            % It writes the combined affine + centering transform to:
            %   <workDir>/full_ct_implant_XFMTO_lpc_mr_pre_do_mat.aff12.1D
            alignScript = fullfile(fileparts(mfilename('fullpath')), ...
                'Utilities', 'eeg_toolbox', 'localize', 'zLocalize', ...
                'shell_scripts', 'align.sh');
            assert(exist(alignScript,'file')==2, ...
                '[electrodeLocalizer] align.sh not found at: %s', alignScript);

            workDir = self.locDirs.ct_1_xfm;
            % Copy inputs to work dir so align.sh can find them by stem
            mrWork = fullfile(workDir, 'mr_pre.nii');
            ctWork = fullfile(workDir, 'ct_implant.nii');
            if exist(mrWork,'file') ~= 2, copyfile(mrNii, mrWork); end
            if exist(ctWork,'file') ~= 2, copyfile(ctNii, ctWork); end

            setenv('PATH', [getenv('PATH') ':' self.afniBin]);
            cmd = sprintf('bash "%s" mr_pre ct_implant "%s"', alignScript, workDir);
            fprintf('[Stage 5] Running: %s\n', cmd);
            [status, txt] = unix(cmd);
            if status ~= 0
                fprintf('%s\n', txt);
                error('[Stage 5] align.sh failed (exit %d).', status);
            end

            % Find the combined transform produced by align.sh.
            % Glob for robustness against naming variations in the script.
            hits = dir(fullfile(workDir, 'full_*.aff12.1D'));
            assert(~isempty(hits), ...
                '[Stage 5] align.sh did not produce a .aff12.1D transform in:\n  %s\nCheck output above.', workDir);
            aff1D = fullfile(workDir, hits(1).name);
            fprintf('[Stage 5] AFNI transform: %s\n', hits(1).name);

            save(xfmFile, 'aff1D');
            fprintf('[Stage 5] Transform path saved to %s\n', xfmFile);
        end

        % -----------------------------------------------------------------
        %% Stage 6 — electrode cluster detection via AFNI
        % -----------------------------------------------------------------

        function detectElectrodes(self, varargin)
            % Detect electrode clusters using AFNI 3dclust (ALICE3 pipeline).
            % Centroids are stored in self.clusters as [N x 3] MR-space mm.
            %
            % Pipeline (mirrors 3dclustering.csh from ALICE3):
            %   1. 3dclust on thresholded CT → labeled cluster volume + clst.1D
            %   2. Resample to 0.5mm + erode/dilate + re-cluster (separates
            %      clusters that touch at CT resolution)
            %   3. Parse clst.1D for centroids in CT RAI mm space
            %   4. xfm_leads to transform centroids to MR space via AFNI

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

            % CT BRIK created by align.sh (3dresample -orient RAI)
            workDir = self.locDirs.ct_1_xfm;
            ctBrik  = fullfile(workDir, 'ct_implant+orig');
            assert(exist([ctBrik '.BRIK'], 'file') == 2, ...
                '[Stage 6] CT BRIK not found: %s.BRIK\n  Re-run coregisterCT.', ctBrik);

            % Ensure AFNI tools are on PATH
            setenv('PATH', [getenv('PATH') ':' self.afniBin]);

            % Clustering subdirectory
            clustDir = fullfile(workDir, 'clustering');
            if ~exist(clustDir, 'dir'), mkdir(clustDir); end

            % Parameters (from 3dclustering.csh defaults)
            rmm  = electrodeLocalizer.AFNI_CLUSTER_RMM;   % neighbor dist in mm
            vmul = electrodeLocalizer.AFNI_CLUSTER_VMUL;  % min volume in µl
            cv   = electrodeLocalizer.CT_HU_THRESHOLD;    % clip threshold (HU)

            outNii = sprintf('3dclusters_r%d_is%d_thr%d.nii', rmm, vmul, cv);
            outFull = fullfile(clustDir, outNii);
            clst1D  = fullfile(clustDir, 'clst.1D');

            fprintf('[Stage 6] Parameters: rmm=%dmm, vmul=%dµl, threshold=%dHU\n', ...
                rmm, vmul, cv);

            % Step 1: initial 3dclust on CT BRIK
            cmd1 = sprintf( ...
                '3dclust -savemask "%s" -overwrite -1Dformat -1clip %d %d %d "%s" > "%s"', ...
                outFull, cv, rmm, vmul, ctBrik, clst1D);
            [status, txt] = unix(sprintf('cd "%s" && %s', clustDir, cmd1));
            if status ~= 0
                fprintf('%s\n', txt);
                error('[Stage 6] 3dclust (step 1) failed (exit %d).', status);
            end

            % Step 2: resample to 0.5mm isotropic, erode -1 then dilate +2,
            %         and re-cluster — separates clusters that touch at CT res.
            tmpRs = 'temp_clusts_rs0.5';
            tmpDe = 'temp_clusts_rs0.5_de2';

            cmd2 = sprintf( ...
                '3dresample -prefix "%s" -overwrite -rmode NN -dxyz 0.5 0.5 0.5 -inset "%s"', ...
                tmpRs, outFull);
            cmd3 = sprintf( ...
                '3dmask_tool -dilate_inputs -1 +2 -prefix "%s" -overwrite -inputs "%s+orig"', ...
                tmpDe, tmpRs);
            cmd4 = sprintf( ...
                '3dclust -savemask "%s" -overwrite -1Dformat -1clip %d %d %d "%s+orig" > "%s"', ...
                outFull, cv, rmm, vmul, tmpDe, clst1D);
            cleanup = sprintf('rm -f "%s+orig.HEAD" "%s+orig.BRIK" "%s+orig.HEAD" "%s+orig.BRIK"', ...
                tmpRs, tmpRs, tmpDe, tmpDe);

            for cmdCell = {cmd2, cmd3, cmd4, cleanup}
                [status, txt] = unix(sprintf('cd "%s" && %s', clustDir, cmdCell{1}));
                if status ~= 0 && ~contains(cmdCell{1}, 'rm')
                    fprintf('%s\n', txt);
                    warning('[Stage 6] AFNI step returned non-zero exit: %s', cmdCell{1}(1:min(60,end)));
                end
            end

            % Step 3: parse clst.1D for cluster centroids (CT RAI mm)
            centroids_ct = electrodeLocalizer.parseClst1D(clst1D);
            N = size(centroids_ct, 1);
            fprintf('[Stage 6] Found %d electrode clusters.\n', N);
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

            self.clusters = clusters_mm;
            save(clustFile, 'clusters_mm');
            fprintf('[Stage 6] %d cluster centroids saved to %s\n', N, clustFile);
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

            dotColors = repmat([0.45 0.45 0.45], N, 1);
            hDots = scatter3(ax, self.clusters(:,1), self.clusters(:,2), ...
                self.clusters(:,3), 50, dotColors, 'filled', 'HitTest', 'off');
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
            xyz      = self.clusters(valid, :);
            self.leads = table(nameOut, xyz(:,1), xyz(:,2), xyz(:,3), typeOut, ...
                'VariableNames', {'chanName','x','y','z','type'});
            fprintf('[Stage 7] %d contacts accepted, %d artifacts.\n', ...
                sum(valid), sum(artifact));

            % ---- nested callbacks (share k, nameOut, etc. by closure) ----

            function refreshDisplay()
                if ~ishandle(fig), return; end
                set(hInfo, 'String', sprintf('Cluster  %d  /  %d', k, N));
                set(hXYZ,  'String', sprintf('x=%.1f   y=%.1f   z=%.1f', ...
                    self.clusters(k,:)));
                set(hHi, 'XData', self.clusters(k,1), ...
                         'YData', self.clusters(k,2), ...
                         'ZData', self.clusters(k,3));

                % List box: remaining (unassigned) channel names
                used = nameOut(~cellfun(@isempty, nameOut));
                if ~isempty(self.chanNames)
                    remaining = self.chanNames(~ismember(self.chanNames, used));
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
                if k > 1, k = k-1; end
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

        function centroids = parseClst1D(clst1DPath)
            % Parse AFNI 3dclust -1Dformat output and return cluster centroids.
            %
            % Returns [N x 3] matrix of centroid coordinates in CT RAI mm
            % (CM_LR, CM_AP, CM_IS — columns 2, 3, 4 of each data row).

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
                return;
            end

            % Each data row: rank  CM_LR  CM_AP  CM_IS  minRL ...
            % Columns 2, 3, 4 (1-indexed) are the RAI mm centroids.
            centroids = zeros(numel(lines), 3);
            for i = 1:numel(lines)
                vals = sscanf(lines{i}, '%f');
                centroids(i, :) = vals(2:4)';
            end
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
