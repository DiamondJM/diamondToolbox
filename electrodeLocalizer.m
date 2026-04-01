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
%   7. manualLocalize    — CT slicer: click to place contacts one-by-one,
%                          then 3-D review (default / recommended); OR
%      namingGUI         — auto-detected cluster naming in 3-D figure
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
        SURFACE_PROXIMITY_MAX_MM = 20;    % max distance (mm) from pial surface to keep a cluster
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
                % Populate chanNames from leads.csv so callers can use them
                if isempty(self.chanNames)
                    leadsFile = fullfile(self.rootFolder, self.subj, 'tal', 'leads.csv');
                    if exist(leadsFile, 'file') == 2
                        T = readtable(leadsFile, 'TextType', 'char');
                        self.chanNames = T.chanName(:);
                    end
                end
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
            dlg = self.localizationSetupDialog(forceNew);
            if strcmp(dlg.action, 'cancel')
                error('electrodeLocalizer:cancelled', ...
                    '[electrodeLocalizer] Setup cancelled by user.');
            end

            % If imports made everything complete (and we're not forcing a
            % re-run), nothing left to do.
            if ~forceNew && self.isComplete(), return; end

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

            % Signal Processing Toolbox
            prereqs.sigproc = license('test','signal_toolbox') && ~isempty(ver('signal'));

            % Statistics and Machine Learning Toolbox
            prereqs.stats = license('test','statistics_toolbox') && ~isempty(ver('stats'));

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

            if ~prereqs.sigproc
                missing{end+1} = sprintf( ...
                    '  [REQUIRED] Signal Processing Toolbox\n%s', ...
                    '             Needed for: EEG time-series processing');
            end

            if ~prereqs.stats
                missing{end+1} = sprintf( ...
                    '  [REQUIRED] Statistics and Machine Learning Toolbox\n%s', ...
                    '             Needed for: clustering and statistical analysis');
            end

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
            electrodeLocalizer.printPrereqLine('Signal Processing Toolbox',           prereqs.sigproc, true);
            electrodeLocalizer.printPrereqLine('Statistics & ML Toolbox',            prereqs.stats,   true);
            electrodeLocalizer.printPrereqLine('gifti toolbox',                      prereqs.gifti,   true);
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
            % CT-to-MR coregistration.
            % Default method is 'auto': AFNI affine alignment.
            % Use 'method','manual' to launch the interactive registration GUI.
            % Saves result to zloc/CT_1/transform/transform.mat.

            p = inputParser;
            addParameter(p, 'forceNew', false);
            addParameter(p, 'method',   'auto');    % 'auto' or 'manual'
            addParameter(p, 'cost', 'lpc');          % auto only: alignment cost function
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;
            method   = p.Results.method;
            cost     = p.Results.cost;

            if strcmp(method, 'manual')
                self.manualCoregisterCT('forceNew', forceNew);
                return;
            end

            xfmFile = fullfile(self.locDirs.ct_1_xfm, 'transform.mat');
            if exist(xfmFile, 'file') == 2 && ~forceNew
                S_chk = load(xfmFile, 'aff1D');
                fi = dir(S_chk.aff1D);
                if ~isempty(fi) && fi.bytes > 0
                    fprintf('[Stage 5] CT-MR transform already exists; skipping.\n');
                    return;
                end
                fprintf('[Stage 5] Stale/empty transform file detected; re-running.\n');
            end

            assert(exist(fullfile(self.afniBin, 'align_epi_anat.py'), 'file') == 2, ...
                '[electrodeLocalizer] AFNI not found at %s. Pass ''afni_bin'' to sourceLocalizer.', self.afniBin);

            fprintf('[Stage 5] Coregistering CT to MR via AFNI align.sh...\n');

            ctNii      = fullfile(self.locDirs.ct_1, 'ct_implant.nii');
            mrNii      = fullfile(self.locDirs.mr_pre, 'mr_pre.nii');
            mrWorkName = 'mr_pre.nii';
            assert(exist(ctNii,'file')==2, '[electrodeLocalizer] CT not found: %s', ctNii);
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

            % Skip align.sh if the combined transform already exists (and is non-empty).
            % align_epi_anat.py can hang on its 3dNotes history step even after
            % the actual alignment finishes; skipping avoids repeated hangs.
            hits = dir(fullfile(workDir, 'full_*.aff12.1D'));
            hits = hits([hits.bytes] > 0);   % ignore stale empty files
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
                unix(sprintf('bash "%s" %s ct_implant "%s" RAI %s > "%s" 2>&1 &', ...
                    alignScript, mrStem, workDir, cost, logFile));

                fprintf('[Stage 5] Running align_epi_anat.py (polling every 15 s) .');
                t0_align = tic; alignTimeout = 3600;
                xfmAf   = []; xfmShft = fullfile(workDir, 'ct_implant_shft.1D');
                while toc(t0_align) < alignTimeout
                    xfmAf = dir(fullfile(workDir, sprintf('*_XFMTO_%s_*_mat.aff12.1D', cost)));
                    % Exclude full_* (created later) and require non-empty file
                    % (3dAllineate touches the file early as a placeholder).
                    xfmAf = xfmAf(~strncmp({xfmAf.name}, 'full_', 5) & [xfmAf.bytes] > 0);
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

            % Enforce rigid-body by stripping any scale/shear from the transform.
            % align_epi_anat.py can introduce scaling even when -rigid_body is
            % requested (due to shell quoting of -Allineate_opts). Polar
            % decomposition via SVD extracts the nearest pure rotation.
            electrodeLocalizer.enforceRigidTransform(aff1D);

            save(xfmFile, 'aff1D');
            fprintf('[Stage 5] Transform path saved to %s\n', xfmFile);
        end


        % -----------------------------------------------------------------
        %% Stage 5 (manual) — FreeView-based manual CT-to-MR coregistration
        % -----------------------------------------------------------------

        function manualCoregisterCT(self, varargin)
            % Launch FreeView for the user to manually align the CT to the MRI.
            % The user should use Tools > Transform Volume in FreeView to align
            % the CT, then File > Save Volume As to save the registered CT as:
            %   <transform_dir>/ct_manual_coreg.nii
            % MATLAB waits for confirmation, then stores the path in transform.mat.

            p = inputParser;
            addParameter(p, 'forceNew', false);
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;

            xfmFile = fullfile(self.locDirs.ct_1_xfm, 'transform.mat');
            if exist(xfmFile, 'file') == 2 && ~forceNew
                S_chk = load(xfmFile);
                if isfield(S_chk, 'manual') && S_chk.manual && ...
                        isfield(S_chk, 'ctCoregNii') && exist(S_chk.ctCoregNii, 'file') == 2
                    fprintf('[Stage 5] Manual CT-MR registration already exists; skipping.\n');
                    return;
                end
            end

            ctNii = fullfile(self.locDirs.ct_1, 'ct_implant.nii');
            mrNii = fullfile(self.locDirs.mr_pre, 'mr_pre.nii');
            assert(exist(ctNii,'file')==2, '[electrodeLocalizer] CT not found: %s', ctNii);
            assert(exist(mrNii,'file')==2, '[electrodeLocalizer] MRI not found: %s', mrNii);

            workDir    = self.locDirs.ct_1_xfm;
            ctWork     = fullfile(workDir, 'ct_implant.nii');
            mrWork     = fullfile(workDir, 'mr_pre.nii');
            ctCoregNii = fullfile(workDir, 'ct_manual_coreg.nii');

            if exist(ctWork, 'file') ~= 2, copyfile(ctNii, ctWork); end
            if exist(mrWork, 'file') ~= 2, copyfile(mrNii, mrWork); end

            fprintf('[Stage 5] Opening CT-MR registration GUI...\n');
            ctCoregNii = self.ctMrRegistrationGUI(mrWork, ctWork, workDir);

            assert(~isempty(ctCoregNii) && exist(ctCoregNii,'file')==2, ...
                '[manualCoregisterCT] Registration cancelled or file not saved.');

            manual = true;  %#ok<NASGU>
            save(xfmFile, 'manual', 'ctCoregNii');
            fprintf('[Stage 5] Manual coregistration saved: %s\n', ctCoregNii);
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

            % Load coregistration result
            workDir = self.locDirs.ct_1_xfm;
            xfmFile = fullfile(workDir, 'transform.mat');
            assert(exist(xfmFile, 'file') == 2, ...
                '[electrodeLocalizer] Run coregisterCT before detectElectrodes.');
            S_xfm = load(xfmFile);

            if isfield(S_xfm, 'manual') && S_xfm.manual
                % Manual coregistration: CT already in MRI space.
                % Resample to RAI for AFNI clustering, use identity transform.
                ctCoregNii = S_xfm.ctCoregNii;
                assert(exist(ctCoregNii,'file')==2, ...
                    '[Stage 6] ct_manual_coreg.nii not found: %s', ctCoregNii);
                ctBrik = fullfile(workDir, 'ct_manual_coreg+orig');
                if ~exist([ctBrik '.BRIK'],'file') && ~exist([ctBrik '.BRIK.gz'],'file')
                    fprintf('[Stage 6] Converting manual coreg CT to AFNI BRIK...\n');
                    setenv('PATH', [getenv('PATH') ':' self.afniBin]);
                    cmd = sprintf('"%s/3dcopy" "%s" "%s"', self.afniBin, ctCoregNii, ctBrik);
                    unix(cmd);
                end
                % Write identity aff12.1D (CT already in MRI space)
                aff1D = fullfile(workDir, 'identity.aff12.1D');
                fid = fopen(aff1D, 'w');
                fprintf(fid, ' 1 0 0 0 0 1 0 0 0 0 1 0\n');
                fclose(fid);
                fprintf('[Stage 6] Using manual coregistration with identity transform.\n');
            else
                % Automatic (AFNI) coregistration
                aff1D = S_xfm.aff1D;
                if exist(aff1D, 'file') ~= 2
                    hits = dir(fullfile(workDir, 'full_*.aff12.1D'));
                    assert(~isempty(hits), ...
                        '[Stage 6] AFNI transform (.aff12.1D) not found. Re-run coregisterCT.');
                    aff1D = fullfile(workDir, hits(1).name);
                    fprintf('[Stage 6] Using transform: %s\n', hits(1).name);
                end
                ctBrik = fullfile(workDir, 'ct_implant+orig');
                assert(exist([ctBrik '.BRIK'], 'file') == 2 || ...
                       exist([ctBrik '.BRIK.gz'], 'file') == 2, ...
                    '[Stage 6] CT BRIK not found: %s.BRIK(.gz)\n  Re-run coregisterCT.', ctBrik);
            end

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
            userQuit   = false;

            % Pre-populate from self.leads if a previous run exists.
            % Fall back to leads.csv on disk if self.leads is empty
            % (e.g. object was re-constructed after pipeline completed).
            if isempty(self.leads)
                leadsFile = fullfile(self.rootFolder, self.subj, 'tal', 'leads.csv');
                if exist(leadsFile, 'file') == 2
                    self.leads = readtable(leadsFile, 'TextType', 'char');
                end
            end
            % Match each cluster to the nearest lead within 15 mm.
            if ~isempty(self.leads)
                leadsXYZ = [self.leads.x, self.leads.y, self.leads.z];
                for ii = 1:N
                    d2 = sum((leadsXYZ - localClusters(ii,:)).^2, 2);
                    [dmin, jj] = min(d2);
                    if sqrt(dmin) < 15
                        nameOut{ii} = char(self.leads.chanName{jj});
                        typeOut{ii} = char(self.leads.type{jj});
                    end
                end
            end

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

            dotColors = repmat([0.25 0.45 0.85], N, 1);   % dark blue — visible against brain surface
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
                'Units','normalized','Position',[x0 0.09 w*0.44 0.07], ...
                'BackgroundColor',[0.30 0.30 0.30],'ForegroundColor',fg, ...
                'FontSize',11,'Callback',@cbBack);

            uicontrol(fig,'Style','pushbutton','String','Finish', ...
                'Units','normalized','Position',[x0+w*0.56 0.09 w*0.44 0.07], ...
                'BackgroundColor',[0.20 0.38 0.58],'ForegroundColor','w', ...
                'FontSize',11,'FontWeight','bold','Callback',@cbFinish);

            uicontrol(fig,'Style','pushbutton','String','Quit', ...
                'Units','normalized','Position',[x0 0.02 w*0.44 0.06], ...
                'BackgroundColor',[0.35 0.15 0.15],'ForegroundColor','w', ...
                'FontSize',11,'Callback',@cbQuit);

            uicontrol(fig,'Style','text', ...
                'String','Drag on brain to rotate', ...
                'Units','normalized','Position',[x0+w*0.56 0.02 w*0.44 0.06], ...
                'BackgroundColor',bg,'ForegroundColor',[0.50 0.50 0.50], ...
                'FontSize',9,'HorizontalAlignment','center');

            refreshDisplay();
            waitfor(fig);   % non-blocking for controls; blocks here until closed

            if userQuit
                error('electrodeLocalizer:cancelled', ...
                    '[electrodeLocalizer] Naming GUI cancelled by user.');
            end

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
                    c(k, :) = [0.25 0.45 0.85];
                    set(hDots, 'CData', c);
                end
                refreshDisplay();
            end

            function cbFinish(~,~)
                delete(fig);
            end

            function cbQuit(~,~)
                choice = questdlg('Quit naming? All assignments will be discarded.', ...
                    'Quit','Quit','Cancel','Cancel');
                if ~strcmp(choice,'Quit'), return; end
                userQuit = true;
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

            workDir     = self.locDirs.ct_1_xfm;
            xfmFilePath = fullfile(workDir, 'transform.mat');
            assert(exist(xfmFilePath,'file')==2, ...
                '[manualLocalize] CT-MR transform not found. Run coregisterCT first.');
            S_xfm = load(xfmFilePath);

            isManualCoreg = isfield(S_xfm, 'manual') && S_xfm.manual;
            if isManualCoreg
                ctFile = S_xfm.ctCoregNii;
                assert(exist(ctFile,'file')==2, ...
                    '[manualLocalize] ct_manual_coreg.nii not found: %s\nRe-run manualCoregisterCT.', ctFile);
            else
                ctFile = fullfile(workDir, 'ct_implant.nii');
                assert(exist(ctFile,'file')==2, ...
                    '[manualLocalize] Coregistered CT not found: %s\nRun coregisterCT first.', ctFile);
                aff1D  = S_xfm.aff1D;
                ctBrik = fullfile(workDir, 'ct_implant+orig');
                assert(exist([ctBrik '.BRIK'],'file')==2 || exist([ctBrik '.BRIK.gz'],'file')==2, ...
                    '[manualLocalize] CT BRIK not found: %s.BRIK(.gz)\nRun coregisterCT first.', ctBrik);
            end

            % ---- Load CT ----
            info = niftiinfo(ctFile);
            vol  = double(niftiread(info));          % [nx, ny, nz]
            Txfm = info.Transform.T;                 % 4×4 MATLAB affine (1-based vox → world)
            [nx, ny, nz] = size(vol);

            % vox2mm: 1-based [i,j,k] row-vector(s) → world mm row-vector(s)
            vox2mm = @(v) ([v, ones(size(v,1),1)] * Txfm);

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

                if isManualCoreg
                    % ---- Manual coreg shortcut ----
                    % CT NIfTI is already in MRI world space (FreeView baked the
                    % user's transform into the header). vox2mm gives FS RAS directly.
                    world_h = vox2mm(vox_0);      % nReal × 4
                    xyzFS   = world_h(:, 1:3);
                else
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

                    % ---- Step 5: MR AFNI DICOM → MR BRIK voxel ----
                    [~, mr_hdr] = unix(sprintf('"%s/3dinfo" -ni -nj -nk -di -dj -dk -o3 "%s" 2>/dev/null', ...
                        self.afniBin, fullfile(workDir,'mr_pre_do+orig')));
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
                    mr_nii_info = niftiinfo(fullfile(workDir, 'mr_pre.nii'));
                    mr_T = mr_nii_info.Transform.T;  % 4×4 MATLAB affine (1-based row-vec → world)
                    mr_vox_1 = mr_nii_vox + 1;       % 0-indexed → 1-indexed
                    mr_ras_h = [mr_vox_1, ones(nReal,1)] * mr_T;  % nReal × 4
                    xyzFS = mr_ras_h(:, 1:3);
                end

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
                lastAx       = [];   % last axes interacted with (for arrow-key scroll)
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
                    'WindowKeyPressFcn',@cbKey, ...
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
                    'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',11,'String',{}, ...
                    'KeyPressFcn',@cbKey);
                uicontrol(fig,'Style','text','String','Or type name:', ...
                    'Units','normalized','Position',[rx 0.635 rw 0.04], ...
                    'BackgroundColor',bg,'ForegroundColor',fg, ...
                    'FontSize',10,'HorizontalAlignment','left');
                hEdit = uicontrol(fig,'Style','edit', ...
                    'Units','normalized','Position',[rx 0.585 rw 0.048], ...
                    'BackgroundColor',[0.25 0.25 0.25],'ForegroundColor',fg, ...
                    'FontSize',11,'HorizontalAlignment','left','String','', ...
                    'KeyPressFcn',@cbKey);
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
                    'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',10,'String',{}, ...
                    'Callback',@cbMkListClick,'KeyPressFcn',@cbKey);
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
                    if ~isempty(ax), lastAx = ax; end
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
                    lastAx   = ax;
                    dragAxes = ax;
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
                            if isempty(ax), ax = lastAx; end
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
                            if isempty(ax), ax = lastAx; end
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

                % ---- Double-click placed marker → jump cursor ----
                function cbMkListClick(~,~)
                    if ~strcmp(get(fig,'SelectionType'),'open'), return; end
                    idx = get(hMkList,'Value');
                    if idx < 1 || idx > numel(markersOut), return; end
                    m = markersOut(idx);
                    invT  = inv(Txfm);
                    vox_h = [m.x, m.y, m.z, 1] * invT;
                    curVox = max(1, min([nx,ny,nz], round(vox_h(1:3))));
                    refreshAll();
                end

                % ---- Remove marker ----
                function cbRemove(~,~)
                    if isempty(markersOut), return; end
                    idx = get(hMkList,'Value');
                    if idx < 1 || idx > numel(markersOut), return; end
                    markersOut(idx) = [];
                    set(hMkList,'Value', max(1, idx-1));
                    pendingPos = curVox;
                    set(hPlaceBtn,'Enable','on');
                    refreshAll();
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

                bg  = [0.12 0.12 0.12];
                fg  = [0.92 0.92 0.92];

                fig3 = figure('Name', sprintf('3-D Review — %s', self.subj), ...
                    'NumberTitle','off','Color',bg, ...
                    'Position',[100 80 1200 820], ...
                    'CloseRequestFcn',@cbQuit);

                ax3 = axes('Parent',fig3,'Position',[0.01 0.01 0.73 0.97], ...
                    'Color','k','XColor','none','YColor','none','ZColor','none');

                mkrXYZ   = [[mkrs.x]', [mkrs.y]', [mkrs.z]'];
                mkrNames = {mkrs.chanName};
                self.renderLeadsOnAxes(ax3, mkrXYZ, mkrNames);

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
        %% Standalone 3-D lead viewer
        % -----------------------------------------------------------------

        function viewLeads(self)
            % Read-only 3-D brain viewer showing named electrode positions.
            % Drag to rotate.  Close the window when done.
            %
            % Usage:  el.viewLeads()

            % Load leads from disk if not in memory.
            if isempty(self.leads)
                leadsFile = fullfile(self.rootFolder, self.subj, 'tal', 'leads.csv');
                assert(exist(leadsFile,'file')==2, ...
                    '[viewLeads] leads.csv not found: %s\nRun the pipeline first.', leadsFile);
                self.leads = readtable(leadsFile, 'TextType', 'char');
            end

            xyz   = [self.leads.x, self.leads.y, self.leads.z];
            names = self.leads.chanName;

            fig = figure('Name', sprintf('Leads — %s', self.subj), ...
                'NumberTitle','off','Color',[0.08 0.08 0.08], ...
                'Position',[80 80 1100 820]);
            ax = axes('Parent',fig,'Color','k', ...
                'XColor','none','YColor','none','ZColor','none', ...
                'Position',[0 0.07 1 0.93]);

            uicontrol('Parent',fig,'Style','pushbutton','String','Close', ...
                'Units','normalized','Position',[0.82 0.01 0.16 0.05], ...
                'BackgroundColor',[0.65 0.10 0.10],'ForegroundColor','w', ...
                'FontSize',11,'Callback',@(~,~)delete(fig));

            self.renderLeadsOnAxes(ax, xyz, names);
            fprintf('[viewLeads] %d electrodes plotted for %s. Close window when done.\n', ...
                size(xyz,1), self.subj);
        end

        % -----------------------------------------------------------------
        %% Shared 3-D brain + electrode renderer
        % -----------------------------------------------------------------

        function renderLeadsOnAxes(self, ax, xyz, names)
            % Render pial brain surfaces + electrode dots + name labels
            % onto an existing axes handle.
            %
            % Inputs:
            %   ax    - axes handle to render into
            %   xyz   - N×3 matrix of electrode positions (FreeSurfer RAS mm)
            %   names - N×1 cell array of electrode names

            sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhFile  = fullfile(sumaDir, 'lh.pial.gii');
            rhFile  = fullfile(sumaDir, 'rh.pial.gii');

            hold(ax,'on'); axis(ax,'equal'); axis(ax,'off');
            view(ax,3); camlight(ax,'headlight'); material(ax,'dull');

            if exist(lhFile,'file')==2
                lhS = gifti(lhFile);
                patch(ax,'Faces',lhS.faces,'Vertices',lhS.vertices, ...
                    'FaceColor',[0.75 0.70 0.65],'EdgeColor','none', ...
                    'FaceAlpha',0.35,'PickableParts','none','HitTest','off');
            end
            if exist(rhFile,'file')==2
                rhS = gifti(rhFile);
                patch(ax,'Faces',rhS.faces,'Vertices',rhS.vertices, ...
                    'FaceColor',[0.75 0.70 0.65],'EdgeColor','none', ...
                    'FaceAlpha',0.35,'PickableParts','none','HitTest','off');
            end

            N = size(xyz,1);
            scatter3(ax, xyz(:,1), xyz(:,2), xyz(:,3), 70, ...
                repmat([0.15 0.35 0.85], N, 1), 'filled', 'HitTest','off', ...
                'MarkerEdgeColor','w','LineWidth',0.5);
            for ii = 1:N
                text(ax, xyz(ii,1)+1, xyz(ii,2), xyz(ii,3), names{ii}, ...
                    'Color',[0.55 0.75 1.00],'FontSize',8, ...
                    'FontWeight','bold','HitTest','off');
            end

            rotate3d(ax,'on');
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

        function ctCoregNii = ctMrRegistrationGUI(~, mrNii, ctNii, outputDir)
        % Interactive rigid-body CT-MR registration GUI.
        % Drag modes: Translate | Rotate | Scale in axial/coronal/sagittal.
        % Save resamples CT onto MRI grid and writes ct_manual_coreg.nii.

        %% ---- Load volumes ----
        mrInfo = niftiinfo(mrNii);
        mrVol  = double(niftiread(mrInfo));
        ctInfo = niftiinfo(ctNii);
        ctVol  = double(niftiread(ctInfo));

        [nxMr, nyMr, nzMr] = size(mrVol);
        pdMr = mrInfo.PixelDimensions(1:3);

        Tmr = mrInfo.Transform.T;
        Tct = ctInfo.Transform.T;

        %% ---- State ----
        tx = 0; ty = 0; tz = 0;  % overwritten below after ctCtrWorld is computed
        rx = 0; ry = 0; rz = 0;
        sc = 1.0;
        ctAlpha = 0.5;
        curVox = [round(nxMr/2), round(nyMr/2), round(nzMr/2)];
        mode_ = 'translate';

        mrWW = 1000; mrWL = 500;
        ctWW = 3000; ctWL = 700;

        ctCtrVox   = [(size(ctVol,1)+1)/2, (size(ctVol,2)+1)/2, (size(ctVol,3)+1)/2];
        ctCtrWorld = [ctCtrVox, 1] * Tct;

        % Initialise translation so CT centre aligns with MRI centre
        mrCtr0 = [round(nxMr/2), round(nyMr/2), round(nzMr/2), 1] * Tmr;
        tx = mrCtr0(1) - ctCtrWorld(1);
        ty = mrCtr0(2) - ctCtrWorld(2);
        tz = mrCtr0(3) - ctCtrWorld(3);

        isDragging   = false;
        dragFigStart = [];
        dragState0   = [];
        dragAxSel    = [];

        SENS_TR  = 0.5;
        SENS_ROT = 0.3;
        SENS_SC  = 0.001;

        %% ---- Orientation flags ----
        colX_ = Tmr(1:3,1);  colY_ = Tmr(1:3,2);  colZ_ = Tmr(1:3,3);
        ax_xFlip  = colX_(1) > 0;
        ax_yFlip  = colY_(2) < 0;
        cor_xFlip = colX_(1) > 0;
        cor_yTop  = colZ_(3) < 0;
        sag_xFlip = colY_(2) < 0;

        %% ---- Colours ----
        bg  = [0.10 0.10 0.10];
        bg2 = [0.18 0.18 0.18];
        bg3 = [0.22 0.22 0.22];
        fg  = [0.92 0.92 0.92];

        %% ---- Figure ----
        ctCoregNii = '';
        fig = figure('Name','CT-MR Registration', ...
            'NumberTitle','off','Color',bg, ...
            'Position',[30 30 1500 870], ...
            'WindowKeyPressFcn',    @cbKey, ...
            'WindowScrollWheelFcn', @cbScroll, ...
            'WindowButtonDownFcn',  @cbDown, ...
            'WindowButtonMotionFcn',@cbMotion, ...
            'WindowButtonUpFcn',    @cbUp, ...
            'CloseRequestFcn',      @cbClose);

        %% ---- Toolbar row 1 ----
        tbH  = 0.048;  tb1Y = 0.950;  tb2Y = tb1Y - tbH - 0.006;

        uicontrol(fig,'Style','text','String','Mode:', ...
            'Units','normalized','Position',[0.01 tb1Y 0.04 tbH], ...
            'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',10);
        btnTr = uicontrol(fig,'Style','togglebutton','String','Translate', ...
            'Units','normalized','Position',[0.055 tb1Y 0.082 tbH], ...
            'BackgroundColor',[0.25 0.55 0.25],'ForegroundColor','w', ...
            'FontSize',10,'Value',1,'Callback',@(~,~)setMode_('translate'));
        btnRot = uicontrol(fig,'Style','togglebutton','String','Rotate', ...
            'Units','normalized','Position',[0.142 tb1Y 0.072 tbH], ...
            'BackgroundColor',bg2,'ForegroundColor',fg, ...
            'FontSize',10,'Value',0,'Callback',@(~,~)setMode_('rotate'));
        btnSc = uicontrol(fig,'Style','togglebutton','String','Scale', ...
            'Units','normalized','Position',[0.219 tb1Y 0.065 tbH], ...
            'BackgroundColor',bg2,'ForegroundColor',fg, ...
            'FontSize',10,'Value',0,'Callback',@(~,~)setMode_('scale'));
        uicontrol(fig,'Style','text', ...
            'String','Drag to transform CT  |  Click to move crosshair  |  Scroll: slice  |  Shift: fine', ...
            'Units','normalized','Position',[0.295 tb1Y 0.44 tbH], ...
            'BackgroundColor',bg,'ForegroundColor',[0.5 0.5 0.5],'FontSize',9, ...
            'HorizontalAlignment','left');

        %% ---- Toolbar row 2 ----
        uicontrol(fig,'Style','text','String','MR W:', ...
            'Units','normalized','Position',[0.01 tb2Y 0.040 tbH], ...
            'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
        hMrW = uicontrol(fig,'Style','edit','String',num2str(mrWW), ...
            'Units','normalized','Position',[0.053 tb2Y 0.055 tbH], ...
            'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',9,'Callback',@cbWL);
        uicontrol(fig,'Style','text','String','L:', ...
            'Units','normalized','Position',[0.110 tb2Y 0.020 tbH], ...
            'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
        hMrL = uicontrol(fig,'Style','edit','String',num2str(mrWL), ...
            'Units','normalized','Position',[0.133 tb2Y 0.055 tbH], ...
            'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',9,'Callback',@cbWL);
        uicontrol(fig,'Style','text','String','CT W:', ...
            'Units','normalized','Position',[0.198 tb2Y 0.040 tbH], ...
            'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
        hCtW = uicontrol(fig,'Style','edit','String',num2str(ctWW), ...
            'Units','normalized','Position',[0.241 tb2Y 0.055 tbH], ...
            'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',9,'Callback',@cbWL);
        uicontrol(fig,'Style','text','String','L:', ...
            'Units','normalized','Position',[0.298 tb2Y 0.020 tbH], ...
            'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
        hCtL = uicontrol(fig,'Style','edit','String',num2str(ctWL), ...
            'Units','normalized','Position',[0.321 tb2Y 0.055 tbH], ...
            'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',9,'Callback',@cbWL);
        uicontrol(fig,'Style','text','String','CT opacity:', ...
            'Units','normalized','Position',[0.387 tb2Y 0.065 tbH], ...
            'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
        hAlphaSlider = uicontrol(fig,'Style','slider','Min',0,'Max',1,'Value',ctAlpha, ...
            'Units','normalized','Position',[0.455 tb2Y+0.005 0.09 tbH*0.8], ...
            'Callback',@cbAlpha);
        addlistener(hAlphaSlider,'Value','PostSet',@(~,~)cbAlphaLive());

        %% ---- Slice panels ----
        figPos  = get(fig,'Position');
        panW    = 0.220;  scrollW = 0.013;
        panH    = panW * figPos(3) / figPos(4);
        panY    = (tb2Y - panH) / 2 + 0.005;

        axAxX  = 0.005;
        slAxX  = axAxX  + panW + 0.001;
        axCorX = slAxX  + scrollW + 0.006;
        slCorX = axCorX + panW + 0.001;
        axSagX = slCorX + scrollW + 0.006;
        slSagX = axSagX + panW + 0.001;

        axAx  = axes('Parent',fig,'Position',[axAxX  panY panW panH], ...
            'Color','k','XColor','none','YColor','none'); hold(axAx,'on');
        axCor = axes('Parent',fig,'Position',[axCorX panY panW panH], ...
            'Color','k','XColor','none','YColor','none'); hold(axCor,'on');
        axSag = axes('Parent',fig,'Position',[axSagX panY panW panH], ...
            'Color','k','XColor','none','YColor','none'); hold(axSag,'on');

        hSlAx  = uicontrol(fig,'Style','slider','Min',1,'Max',nzMr,'Value',curVox(3), ...
            'SliderStep',[1/(nzMr-1), 10/(nzMr-1)], ...
            'Units','normalized','Position',[slAxX  panY scrollW panH], ...
            'BackgroundColor',bg2,'Callback',@(~,~)cbSlider(1));
        hSlCor = uicontrol(fig,'Style','slider','Min',1,'Max',nyMr,'Value',curVox(2), ...
            'SliderStep',[1/(nyMr-1), 10/(nyMr-1)], ...
            'Units','normalized','Position',[slCorX panY scrollW panH], ...
            'BackgroundColor',bg2,'Callback',@(~,~)cbSlider(2));
        hSlSag = uicontrol(fig,'Style','slider','Min',1,'Max',nxMr,'Value',curVox(1), ...
            'SliderStep',[1/(nxMr-1), 10/(nxMr-1)], ...
            'Units','normalized','Position',[slSagX panY scrollW panH], ...
            'BackgroundColor',bg2,'Callback',@(~,~)cbSlider(3));
        addlistener(hSlAx,  'Value','PostSet',@(~,~)cbSlider(1));
        addlistener(hSlCor, 'Value','PostSet',@(~,~)cbSlider(2));
        addlistener(hSlSag, 'Value','PostSet',@(~,~)cbSlider(3));

        hImAx  = []; hImCor = []; hImSag = [];
        hXHax  = []; hXHcor = []; hXHsag = [];
        hTitleAx  = title(axAx,  '','Color',fg,'FontSize',9);
        hTitleCor = title(axCor, '','Color',fg,'FontSize',9);
        hTitleSag = title(axSag, '','Color',fg,'FontSize',9);

        %% ---- Status bar ----
        hStatus = uicontrol(fig,'Style','text', ...
            'Units','normalized','Position',[0.005 0.003 0.73 0.038], ...
            'BackgroundColor',bg,'ForegroundColor',[0.5 0.8 1.0], ...
            'FontSize',9,'HorizontalAlignment','left','String','');

        %% ---- Right panel ----
        rp = 0.742;  rw = 0.253;
        uicontrol(fig,'Style','text','String','CT-MR Registration', ...
            'Units','normalized','Position',[rp 0.900 rw 0.048], ...
            'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',12,'FontWeight','bold');

        labels_  = {'Tx (mm):','Ty (mm):','Tz (mm):','Rx (°):','Ry (°):','Rz (°):','Scale:'};
        hReadout = gobjects(1,7);
        for ii = 1:7
            yy = 0.840 - (ii-1)*0.052;
            uicontrol(fig,'Style','text','String',labels_{ii}, ...
                'Units','normalized','Position',[rp yy 0.095 0.038], ...
                'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',10, ...
                'HorizontalAlignment','right');
            hReadout(ii) = uicontrol(fig,'Style','edit','String','0.0', ...
                'Units','normalized','Position',[rp+0.100 yy 0.080 0.038], ...
                'BackgroundColor',bg2,'ForegroundColor',[0.3 0.9 0.3],'FontSize',10, ...
                'HorizontalAlignment','center', ...
                'Callback',@(h,~) cbEditParam(ii, str2double(get(h,'String'))));
        end

        uicontrol(fig,'Style','pushbutton','String','Reset to Center', ...
            'Units','normalized','Position',[rp 0.208 rw*0.95 0.048], ...
            'BackgroundColor',[0.55 0.25 0.25],'ForegroundColor','w', ...
            'FontSize',10,'Callback',@cbReset);
        uicontrol(fig,'Style','pushbutton','String','Save', ...
            'Units','normalized','Position',[rp 0.100 rw*0.44 0.065], ...
            'BackgroundColor',[0.2 0.5 0.2],'ForegroundColor','w', ...
            'FontSize',12,'FontWeight','bold','Callback',@cbSave);
        uicontrol(fig,'Style','pushbutton','String','Cancel', ...
            'Units','normalized','Position',[rp+rw*0.51 0.100 rw*0.44 0.065], ...
            'BackgroundColor',[0.5 0.25 0.25],'ForegroundColor','w', ...
            'FontSize',12,'Callback',@cbClose);

        %% ---- Initial draw ----
        refreshAll();
        uiwait(fig);

        % ==============================================================
        % Nested functions
        % ==============================================================

            function T = buildRegTransform()
                ctr = ctCtrWorld(1:3);
                Tc  = eye(4);  Tc(4,1:3)  = -ctr;
                Tci = eye(4);  Tci(4,1:3) =  ctr;
                Ttr = eye(4);  Ttr(4,1:3) = [tx ty tz];
                S   = diag([sc sc sc 1]);
                cx = cosd(rx); sx_ = sind(rx);
                cy = cosd(ry); sy_ = sind(ry);
                cz = cosd(rz); sz_ = sind(rz);
                Rx = [1  0    0   0; 0  cx   sx_ 0; 0 -sx_  cx  0; 0  0    0   1];
                Ry = [cy  0  -sy_ 0; 0   1   0   0; sy_ 0   cy  0; 0   0   0   1];
                Rz = [cz  sz_ 0  0; -sz_ cz  0  0;  0   0   1  0;  0   0   0   1];
                T = Tc * S * Rx * Ry * Rz * Tci * Ttr;
            end

            function M = buildSamplingMatrix()
                T_reg = buildRegTransform();
                M = Tmr / T_reg / Tct;
            end

            function ctSlice = sampleCTslice(dim, sliceIdx)
                M = buildSamplingMatrix();
                switch dim
                    case 1
                        [jj, ii] = meshgrid(1:nyMr, 1:nxMr);
                        kk = repmat(sliceIdx, nxMr, nyMr);
                    case 2
                        [kk, ii] = meshgrid(1:nzMr, 1:nxMr);
                        jj = repmat(sliceIdx, nxMr, nzMr);
                    case 3
                        [kk, jj] = meshgrid(1:nzMr, 1:nyMr);
                        ii = repmat(sliceIdx, nyMr, nzMr);
                end
                sz_out = size(ii);  N = numel(ii);
                mr_h = [ii(:), jj(:), kk(:), ones(N,1)];
                ct_h = mr_h * M;
                ctSlice = reshape(interp3(ctVol, ct_h(:,2), ct_h(:,1), ct_h(:,3), ...
                    'linear', NaN), sz_out);
            end

            function rgb = makeOverlay(mrSlice, ctSlice)
                mr_lo = mrWL - mrWW/2;  mr_hi = mrWL + mrWW/2;
                ct_lo = ctWL - ctWW/2;  ct_hi = ctWL + ctWW/2;
                mr_n  = min(max((mrSlice - mr_lo)/(mr_hi - mr_lo), 0), 1);
                mr_rgb = repmat(mr_n, 1, 1, 3);
                if isempty(ctSlice) || all(isnan(ctSlice(:))), rgb = mr_rgb; return; end
                ct_n  = min(max((ctSlice - ct_lo)/(ct_hi - ct_lo), 0), 1);
                hmap  = hot(256);
                idx   = max(1, min(256, floor(ct_n*255)+1));
                ct_rgb = cat(3, reshape(hmap(idx(:),1),size(ct_n)), ...
                                reshape(hmap(idx(:),2),size(ct_n)), ...
                                reshape(hmap(idx(:),3),size(ct_n)));
                a  = ctAlpha * sqrt(ct_n) .* double(~isnan(ctSlice));
                a3 = repmat(a, 1, 1, 3);
                rgb = min(max((1-a3).*mr_rgb + a3.*ct_rgb, 0), 1);
            end

            function refreshAll()
                if ~ishandle(fig), return; end
                % Axial
                sl_mr  = mrVol(:,:,curVox(3))';
                ct_raw = sampleCTslice(1, curVox(3))';
                if ax_yFlip, sl_mr = flipud(sl_mr); ct_raw = flipud(ct_raw); end
                if isempty(hImAx) || ~ishandle(hImAx)
                    hImAx = image(axAx, makeOverlay(sl_mr, ct_raw));
                    axis(axAx,'normal');
                    set(axAx,'XLim',[0.5 nxMr+0.5],'YLim',[0.5 nyMr+0.5],'YDir','normal', ...
                        'DataAspectRatio',[pdMr(1) pdMr(2) 1]);
                    if ax_xFlip, set(axAx,'XDir','reverse'); end
                    addOrientLabels(axAx,'R','L','A','P');
                else
                    set(hImAx,'CData', makeOverlay(sl_mr, ct_raw));
                end
                set(hTitleAx,'String',sprintf('AXIAL  z=%d/%d',curVox(3),nzMr));
                % Coronal
                sl_mr  = squeeze(mrVol(:,curVox(2),:))';
                ct_raw = sampleCTslice(2, curVox(2))';
                if cor_yTop, sl_mr = flipud(sl_mr); ct_raw = flipud(ct_raw); end
                if isempty(hImCor) || ~ishandle(hImCor)
                    hImCor = image(axCor, makeOverlay(sl_mr, ct_raw));
                    axis(axCor,'normal');
                    set(axCor,'XLim',[0.5 nxMr+0.5],'YLim',[0.5 nzMr+0.5],'YDir','normal', ...
                        'DataAspectRatio',[pdMr(1) pdMr(3) 1]);
                    if cor_xFlip, set(axCor,'XDir','reverse'); end
                    addOrientLabels(axCor,'R','L','S','I');
                else
                    set(hImCor,'CData', makeOverlay(sl_mr, ct_raw));
                end
                set(hTitleCor,'String',sprintf('CORONAL  y=%d/%d',curVox(2),nyMr));
                % Sagittal
                sl_mr  = squeeze(mrVol(curVox(1),:,:))';
                ct_raw = sampleCTslice(3, curVox(1))';
                if cor_yTop, sl_mr = flipud(sl_mr); ct_raw = flipud(ct_raw); end
                if isempty(hImSag) || ~ishandle(hImSag)
                    hImSag = image(axSag, makeOverlay(sl_mr, ct_raw));
                    axis(axSag,'normal');
                    set(axSag,'XLim',[0.5 nyMr+0.5],'YLim',[0.5 nzMr+0.5],'YDir','normal', ...
                        'DataAspectRatio',[pdMr(2) pdMr(3) 1]);
                    if sag_xFlip, set(axSag,'XDir','reverse'); end
                    addOrientLabels(axSag,'P','A','S','I');
                else
                    set(hImSag,'CData', makeOverlay(sl_mr, ct_raw));
                end
                set(hTitleSag,'String',sprintf('SAGITTAL  x=%d/%d',curVox(1),nxMr));
                refreshCrosshairs();
                updateReadouts();
                if ishandle(hSlAx),  set(hSlAx,  'Value', curVox(3)); end
                if ishandle(hSlCor), set(hSlCor, 'Value', curVox(2)); end
                if ishandle(hSlSag), set(hSlSag, 'Value', curVox(1)); end
            end

            function refreshCrosshairs()
                if ~ishandle(fig), return; end
                xhColor = [0.3 1.0 0.3];  lw = 0.8;
                iy_ax  = electrodeLocalizer.ternary(ax_yFlip,  nyMr-curVox(2)+1, curVox(2));
                iz_cor = electrodeLocalizer.ternary(cor_yTop, nzMr-curVox(3)+1, curVox(3));
                if ~isempty(hXHax)  && all(ishandle(hXHax)),  delete(hXHax);  end
                if ~isempty(hXHcor) && all(ishandle(hXHcor)), delete(hXHcor); end
                if ~isempty(hXHsag) && all(ishandle(hXHsag)), delete(hXHsag); end
                hXHax(1)  = plot(axAx, [curVox(1) curVox(1)],[0.5 nyMr+0.5], 'Color',xhColor,'LineWidth',lw,'HitTest','off');
                hXHax(2)  = plot(axAx, [0.5 nxMr+0.5],[iy_ax iy_ax],         'Color',xhColor,'LineWidth',lw,'HitTest','off');
                hXHcor(1) = plot(axCor,[curVox(1) curVox(1)],[0.5 nzMr+0.5], 'Color',xhColor,'LineWidth',lw,'HitTest','off');
                hXHcor(2) = plot(axCor,[0.5 nxMr+0.5],[iz_cor iz_cor],        'Color',xhColor,'LineWidth',lw,'HitTest','off');
                hXHsag(1) = plot(axSag,[curVox(2) curVox(2)],[0.5 nzMr+0.5], 'Color',xhColor,'LineWidth',lw,'HitTest','off');
                hXHsag(2) = plot(axSag,[0.5 nyMr+0.5],[iz_cor iz_cor],        'Color',xhColor,'LineWidth',lw,'HitTest','off');
            end

            function updateReadouts()
                vals = {tx,ty,tz,rx,ry,rz,sc};
                fmts = {'%.1f','%.1f','%.1f','%.1f','%.1f','%.1f','%.4f'};
                for jj = 1:7
                    set(hReadout(jj),'String',sprintf(fmts{jj}, vals{jj}));
                end
                set(hStatus,'String',sprintf( ...
                    '  T=[%.1f, %.1f, %.1f] mm   R=[%.1f, %.1f, %.1f]°   scale=%.4f', ...
                    tx,ty,tz,rx,ry,rz,sc));
            end

            function cbDown(~,~)
                ax = axUnderCursor();
                if isempty(ax), return; end
                isDragging = false;  dragFigStart = get(fig,'CurrentPoint');
                dragAxSel  = ax;     dragState0   = [tx ty tz rx ry rz sc];
            end

            function cbMotion(~,~)
                if isempty(dragFigStart), return; end
                delta = get(fig,'CurrentPoint') - dragFigStart;
                if ~isDragging && norm(delta) < 3, return; end
                isDragging = true;
                fine = any(strcmp(get(fig,'CurrentModifier'),'shift'));
                k    = electrodeLocalizer.ternary(fine, 0.1, 1.0);
                dx   = delta(1);  dy = delta(2);
                xSignAx = -1;  xSignSag = 1;
                ax = dragAxSel;
                switch mode_
                    case 'translate'
                        s = k*SENS_TR;
                        if     ax==axAx,  tx=dragState0(1)+dx*s*xSignAx;  ty=dragState0(2)+dy*s;
                        elseif ax==axCor, tx=dragState0(1)+dx*s*xSignAx;  tz=dragState0(3)+dy*s;
                        elseif ax==axSag, ty=dragState0(2)+dx*s*xSignSag; tz=dragState0(3)+dy*s;
                        end
                    case 'rotate'
                        s = k*SENS_ROT;
                        if     ax==axAx,  rz=dragState0(6)+dx*s*(-xSignAx);
                        elseif ax==axCor, ry=dragState0(5)+dx*s*xSignAx;
                        elseif ax==axSag, rx=dragState0(4)+dx*s*(-xSignSag);
                        end
                    case 'scale'
                        sc = max(0.1, dragState0(7)+dx*k*SENS_SC);
                end
                refreshAll();
            end

            function cbUp(~,~)
                if ~isDragging && ~isempty(dragFigStart) && ~isempty(dragAxSel)
                    ax = dragAxSel;
                    cp = get(ax,'CurrentPoint');
                    xd = round(cp(1,1));  yd = round(cp(1,2));
                    if ax==axAx
                        curVox(1) = max(1,min(nxMr,xd));
                        curVox(2) = electrodeLocalizer.ternary(ax_yFlip, nyMr-max(1,min(nyMr,yd))+1, max(1,min(nyMr,yd)));
                    elseif ax==axCor
                        curVox(1) = max(1,min(nxMr,xd));
                        curVox(3) = electrodeLocalizer.ternary(cor_yTop, nzMr-max(1,min(nzMr,yd))+1, max(1,min(nzMr,yd)));
                    elseif ax==axSag
                        curVox(2) = max(1,min(nyMr,xd));
                        curVox(3) = electrodeLocalizer.ternary(cor_yTop, nzMr-max(1,min(nzMr,yd))+1, max(1,min(nzMr,yd)));
                    end
                    refreshAll();
                end
                isDragging=false; dragFigStart=[]; dragAxSel=[]; dragState0=[];
            end

            function cbScroll(~,evt)
                ax = axUnderCursor();  d = evt.VerticalScrollCount;
                if     ax==axAx,  curVox(3)=max(1,min(nzMr,curVox(3)-d));
                elseif ax==axCor, curVox(2)=max(1,min(nyMr,curVox(2)-d));
                elseif ax==axSag, curVox(1)=max(1,min(nxMr,curVox(1)-d));
                end
                refreshAll();
            end

            function cbKey(~,evt)
                switch evt.Key
                    case 'uparrow',   curVox(3)=min(nzMr,curVox(3)+1); refreshAll();
                    case 'downarrow', curVox(3)=max(1,   curVox(3)-1); refreshAll();
                end
            end

            function setMode_(m)
                mode_ = m;
                set(btnTr,  'Value',strcmp(m,'translate'), 'BackgroundColor', ...
                    electrodeLocalizer.ternary(strcmp(m,'translate'),[0.25 0.55 0.25],bg2));
                set(btnRot, 'Value',strcmp(m,'rotate'),    'BackgroundColor', ...
                    electrodeLocalizer.ternary(strcmp(m,'rotate'),   [0.25 0.55 0.25],bg2));
                set(btnSc,  'Value',strcmp(m,'scale'),     'BackgroundColor', ...
                    electrodeLocalizer.ternary(strcmp(m,'scale'),    [0.25 0.55 0.25],bg2));
            end

            function cbWL(~,~)
                mrWW=str2double(get(hMrW,'String')); mrWL=str2double(get(hMrL,'String'));
                ctWW=str2double(get(hCtW,'String')); ctWL=str2double(get(hCtL,'String'));
                refreshAll();
            end

            function cbAlpha(~,~),     ctAlpha=get(hAlphaSlider,'Value'); refreshAll(); end
            function cbAlphaLive(),    ctAlpha=get(hAlphaSlider,'Value'); refreshAll(); end

            function cbEditParam(idx, val)
                if isnan(val), return; end
                switch idx
                    case 1,tx=val; case 2,ty=val; case 3,tz=val;
                    case 4,rx=val; case 5,ry=val; case 6,rz=val;
                    case 7,sc=max(0.1,val);
                end
                refreshAll();
            end

            function cbSlider(dim)
                switch dim
                    case 1, curVox(3)=round(get(hSlAx, 'Value'));
                    case 2, curVox(2)=round(get(hSlCor,'Value'));
                    case 3, curVox(1)=round(get(hSlSag,'Value'));
                end
                refreshAll();
            end

            function cbReset(~,~)
                mrCtr = [round(nxMr/2),round(nyMr/2),round(nzMr/2),1]*Tmr;
                tx=mrCtr(1)-ctCtrWorld(1); ty=mrCtr(2)-ctCtrWorld(2); tz=mrCtr(3)-ctCtrWorld(3);
                rx=0; ry=0; rz=0; sc=1;
                refreshAll();
            end

            function cbSave(~,~)
                set(hStatus,'String','  Resampling CT onto MRI grid...'); drawnow;
                M  = buildSamplingMatrix();
                [ii,jj,kk] = ndgrid(1:nxMr,1:nyMr,1:nzMr);  N = numel(ii);
                ct_h = [ii(:),jj(:),kk(:),ones(N,1)] * M;
                ctRs = reshape(interp3(ctVol,ct_h(:,2),ct_h(:,1),ct_h(:,3),'linear',0), ...
                    nxMr,nyMr,nzMr);
                outPath = fullfile(outputDir,'ct_manual_coreg.nii');
                outInfo = mrInfo;  outInfo.Datatype='single';  outInfo.BitsPerPixel=32;
                niftiwrite(single(ctRs), outPath, outInfo, 'Compressed', false);
                ctCoregNii = outPath;
                fprintf('[ctMrRegistrationGUI] Saved: %s\n', outPath);
                uiresume(fig); delete(fig);
            end

            function cbClose(~,~)
                ctCoregNii = '';
                if ishandle(fig), uiresume(fig); delete(fig); end
            end

            function ax = axUnderCursor()
                ax = [];
                fp  = get(fig,'CurrentPoint');  fsz = get(fig,'Position');
                np  = fp ./ fsz(3:4);
                for aa = [axAx axCor axSag]
                    apos = get(aa,'Position');
                    if np(1)>=apos(1) && np(1)<=apos(1)+apos(3) && ...
                       np(2)>=apos(2) && np(2)<=apos(2)+apos(4)
                        ax = aa; return;
                    end
                end
            end

            function addOrientLabels(ax, xl, xr, yt, yb)
                xl_ = get(ax,'XLim');  yl_ = get(ax,'YLim');
                mx = mean(xl_);  my = mean(yl_);
                text(ax,xl_(1)+0.01*(xl_(2)-xl_(1)),my,xl,'Color',[1 0.8 0],'FontSize',8, ...
                    'HorizontalAlignment','left','VerticalAlignment','middle','HitTest','off');
                text(ax,xl_(2)-0.01*(xl_(2)-xl_(1)),my,xr,'Color',[1 0.8 0],'FontSize',8, ...
                    'HorizontalAlignment','right','VerticalAlignment','middle','HitTest','off');
                text(ax,mx,yl_(2)-0.01*(yl_(2)-yl_(1)),yt,'Color',[1 0.8 0],'FontSize',8, ...
                    'HorizontalAlignment','center','VerticalAlignment','top','HitTest','off');
                text(ax,mx,yl_(1)+0.01*(yl_(2)-yl_(1)),yb,'Color',[1 0.8 0],'FontSize',8, ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom','HitTest','off');
            end

        end % ctMrRegistrationGUI

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

        function dlg = localizationSetupDialog(self, forceNew)
            % Dark-themed modal dialog shown when required localization files
            % are missing, or when forceNew=true.  A listbox shows all 5
            % files with OK/blank status; selecting a row updates the
            % description panel.  Create runs the full pipeline; Import
            % copies an existing file into place (always allowed, even if
            % the file is already present).
            %
            % Returns struct with .action: 'create' | 'import' | 'cancel'
            if nargin < 2, forceNew = false; end

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
            if forceNew
                hdrTitle = sprintf('Re-run localization for  %s', self.subj);
                hdrSub   = 'forceNew=true: Create re-runs the full pipeline.  Import replaces any file.';
            else
                hdrTitle = sprintf('Localization files not found for  %s', self.subj);
                hdrSub   = 'Select a row to see details.  Create runs the full pipeline.  Import copies existing files.';
            end
            uicontrol(fig, 'Style','text', ...
                'String', hdrTitle, ...
                'ForegroundColor', FG, 'BackgroundColor', BG, ...
                'FontSize', 12, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'left', ...
                'Position', [PAD hdrY+28 W-PAD*2 22]);
            uicontrol(fig, 'Style','text', ...
                'String', hdrSub, ...
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
                % Import the currently selected file; update its row
                % in-place and stay on the dialog.  Re-import is always
                % allowed so forceNew users can replace existing files.
                k = get(hList, 'Value');
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

        function enforceRigidTransform(aff1DPath)
            % Read an AFNI aff12.1D, strip any scale/shear via SVD polar
            % decomposition, and write back a pure rotation + translation.
            fid = fopen(aff1DPath, 'r');
            C   = textscan(fid, '%f');
            fclose(fid);
            T = C{1}(:)';                    % 1×12 row vector
            assert(numel(T) == 12, ...
                'enforceRigidTransform: expected 12 values in %s, got %d', aff1DPath, numel(T));
            % aff12.1D layout: r11 r12 r13 dx  r21 r22 r23 dy  r31 r32 r33 dz
            A = [T(1:3); T(5:7); T(9:11)];  % 3×3 rotation matrix (one row per line)
            t = [T(4), T(8), T(12)];         % 1×3 translation
            [U, ~, V] = svd(A);
            R = U * V';
            if det(R) < 0                    % avoid improper rotation (reflection)
                V(:,3) = -V(:,3);
                R = U * V';
            end
            T_rigid = [R(1,:), t(1), R(2,:), t(2), R(3,:), t(3)];
            fid = fopen(aff1DPath, 'w');
            fprintf(fid, ' %12.8f', T_rigid);
            fprintf(fid, '\n');
            fclose(fid);
            fprintf('[Stage 5] Rigid-body transform enforced (scale stripped).\n');
        end

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
