
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
%   5. coregisterCT      — SPM rigid-body CT→MR registration
%   6. detectElectrodes  — threshold + connected-components on CT volume
%   7. namingGUI         — interactive 3-D figure; user names each cluster,
%                          labels it depth or subdural, or marks as artifact
%   8. projectElectrodes — subdural contacts snapped to nearest pial vertex;
%                          depth contacts kept at CT-MR coordinates
%   9. writeLeads        — write tal/leads.csv (chanName, x, y, z)
%
% DEPENDENCIES
%   FreeSurfer  — recon-all surface reconstruction
%   AFNI/SUMA   — @SUMA_Make_Spec_FS standard mesh (ld141, 198 812 vertices)
%   SPM         — spm_coreg CT-to-MR coregistration (must be on MATLAB path)
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
        CT_HU_THRESHOLD  = 1500;   % HU threshold for electrode detection
        MIN_CLUSTER_VOXELS = 3;    % minimum voxels per cluster
        MAX_CLUSTER_VOXELS = 500;  % maximum voxels per cluster (reject large artifacts)
        AFNI_BIN_DEFAULT = fullfile(char(java.lang.System.getProperty('user.home')), 'abin');
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

            self.run();
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
                    sprintf('Localization files (leads.csv, lh.pial.gii, rh.pial.gii) not found for %s.\n\nHow would you like to proceed?', self.subj), ...
                    'electrodeLocalizer', ...
                    'Run localization pipeline to create them', 'These files already exist; import them', 'Cancel', ...
                    'Run localization pipeline to create them');

                if isempty(choice) || strcmp(choice, 'Cancel')
                    fprintf('[electrodeLocalizer] Cancelled.\n');
                    return;
                end

                if strcmp(choice, 'Import Existing Files')
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
            self.runSurface('forceNew', forceNew);
            self.runSuma('forceNew', forceNew);
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
            %   .spm
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
            xfmFile    = fullfile(self.locDirs.ct_1_xfm, 'transform_spm.mat');
            clustFile  = fullfile(self.locDirs.ct_1_xfm, 'clusters_mr.mat');

            needSurface = ~(exist(lhPial,'file')==2 && exist(rhPial,'file')==2);
            needSuma    = ~(exist(lhGii,'file')==2  && exist(rhGii,'file')==2);
            needCoreg   = exist(xfmFile, 'file') ~= 2;
            needDetect  = exist(clustFile,'file') ~= 2;

            % --- Check each prerequisite ---

            % Image Processing Toolbox (bwconncomp, regionprops3)
            prereqs.imageProcessingToolbox = ...
                license('test','image_toolbox') && exist('bwconncomp','file') == 2;

            % SPM
            prereqs.spm = exist('spm_coreg','file') == 2;

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

            if ~prereqs.imageProcessingToolbox && needDetect
                missing{end+1} = sprintf( ...
                    '  [REQUIRED] MATLAB Image Processing Toolbox\n%s\n%s', ...
                    '             Needed for: electrode cluster detection (stage 6)', ...
                    '             Fix: install the Image Processing Toolbox add-on.');
            end

            if ~prereqs.spm && needCoreg
                missing{end+1} = sprintf( ...
                    '  [REQUIRED] SPM\n%s\n%s', ...
                    '             Needed for: CT-MR coregistration (stage 5)', ...
                    '             Fix: download SPM (https://www.fil.ion.ucl.ac.uk/spm/) and add to MATLAB path.');
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
            electrodeLocalizer.printPrereqLine('Image Processing Toolbox', prereqs.imageProcessingToolbox, needDetect);
            electrodeLocalizer.printPrereqLine('SPM (spm_coreg)',           prereqs.spm,                   needCoreg);
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

            surfDir  = fullfile(self.locDirs.fs_subj, 'surf');
            lhPial   = fullfile(surfDir, 'lh.pial-outer-smoothed');
            rhPial   = fullfile(surfDir, 'rh.pial-outer-smoothed');
            origMgz  = fullfile(self.locDirs.fs_subj, 'mri', 'orig.mgz');
            done    = exist(lhPial, 'file') == 2 && exist(rhPial, 'file') == 2 && exist(origMgz, 'file') == 2;

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
            suma(self.subj, self.locDirs.fs, ...
                'afni_bin',       self.afniBin, ...
                'freesurfer_bin', self.fsBin, ...
                'rerun',          forceNew);
        end

        % -----------------------------------------------------------------
        %% Stage 5 — CT-to-MR coregistration via SPM
        % -----------------------------------------------------------------

        function coregisterCT(self, varargin)
            % Rigid-body coregistration of post-op CT to pre-op MRI using
            % SPM's normalised mutual information cost function.
            % Saves the resulting 4x4 world-space transform to
            % zloc/CT_1/transform/transform_spm.mat.

            p = inputParser;
            addParameter(p, 'forceNew', false);
            parse(p, varargin{:});
            forceNew = p.Results.forceNew;

            xfmFile = fullfile(self.locDirs.ct_1_xfm, 'transform_spm.mat');
            if exist(xfmFile, 'file') == 2 && ~forceNew
                fprintf('[Stage 5] CT-MR transform already exists; skipping.\n');
                return;
            end

            assert(exist('spm_coreg', 'file') == 2, ...
                '[electrodeLocalizer] SPM not found on MATLAB path. Please add SPM.');

            fprintf('[Stage 5] Coregistering CT to MR via SPM...\n');

            mrNii = fullfile(self.locDirs.mr_pre,  'mr_pre.nii');
            ctNii = fullfile(self.locDirs.ct_1,    'ct_implant.nii');

            VF = spm_vol(mrNii);   % fixed  = MRI
            VG = spm_vol(ctNii);   % moving = CT

            % spm_coreg returns 6 rigid-body params [x y z pitch roll yaw]
            % in the world coordinate frames of both volumes.
            flags.cost_fun = 'nmi';
            flags.sep      = [4 2];
            flags.tol      = [0.02 0.02 0.02 0.001 0.001 0.001];
            flags.fwhm     = [7 7];
            x = spm_coreg(VF, VG, flags);

            % Build the 4x4 transform: CT voxels → MRI world mm (RAS).
            %
            % spm_coreg(VF, VG) returns x such that spm_matrix(x) maps
            % CT world coordinates into MRI world coordinates.  Prepending
            % VG.mat (CT vox-to-world) gives a single matrix that takes
            % CT voxel indices directly to MRI RAS mm:
            %
            %   P_mr_mm = M_ct2mr * [vox_i; vox_j; vox_k; 1]
            %
            M_ct2mr = spm_matrix(x) * VG.mat;

            save(xfmFile, 'M_ct2mr', 'x');
            fprintf('[Stage 5] Transform saved to %s\n', xfmFile);
        end

        % -----------------------------------------------------------------
        %% Stage 6 — electrode cluster detection
        % -----------------------------------------------------------------

        function detectElectrodes(self, varargin)
            % Threshold the CT volume and extract cluster centroids.
            % Centroids are stored in self.clusters as [N x 3] MR-space mm.

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

            fprintf('[Stage 6] Detecting electrode clusters in CT...\n');

            % Load coregistration transform
            xfmFile = fullfile(self.locDirs.ct_1_xfm, 'transform_spm.mat');
            assert(exist(xfmFile, 'file') == 2, ...
                '[electrodeLocalizer] Run coregisterCT before detectElectrodes.');
            S = load(xfmFile, 'M_ct2mr');

            % Load CT volume
            ctNii = fullfile(self.locDirs.ct_1, 'ct_implant.nii');
            V  = spm_vol(ctNii);
            ct = spm_read_vols(V);

            % Threshold and connected-component analysis
            mask = ct >= electrodeLocalizer.CT_HU_THRESHOLD;
            CC   = bwconncomp(mask, 26);   % 26-connectivity in 3D
            stats = regionprops3(CC, 'Centroid', 'Volume');

            % Filter by voxel count
            keep = stats.Volume >= electrodeLocalizer.MIN_CLUSTER_VOXELS & ...
                   stats.Volume <= electrodeLocalizer.MAX_CLUSTER_VOXELS;
            centroids_vox = stats.Centroid(keep, :);   % [N x 3] voxel indices

            fprintf('[Stage 6] Found %d candidate electrode clusters.\n', sum(keep));

            % Convert CT voxel centroids directly to MRI RAS mm.
            % M_ct2mr = spm_matrix(x) * VG.mat already encodes the full
            % CT-vox → MRI-world chain, so no intermediate step needed.
            N   = size(centroids_vox, 1);
            hom = [centroids_vox'; ones(1, N)];   % 4 x N homogeneous
            mr_world    = S.M_ct2mr * hom;        % MRI RAS mm (4 x N)
            clusters_mm = mr_world(1:3, :)';      % N x 3

            self.clusters = clusters_mm;
            save(clustFile, 'clusters_mm');
            fprintf('[Stage 6] Cluster centroids saved to %s\n', clustFile);
        end

        % -----------------------------------------------------------------
        %% Stage 7 — interactive naming GUI
        % -----------------------------------------------------------------

        function namingGUI(self)
            % Display detected clusters on the pial surface and let the
            % user name each contact, specify depth vs subdural, or mark
            % it as artifact.
            %
            % Results stored in self.leads as a table with columns:
            %   chanName, x, y, z, type ('depth' | 'subdural')

            assert(~isempty(self.clusters), ...
                '[electrodeLocalizer] Run detectElectrodes before namingGUI.');

            N = size(self.clusters, 1);
            fprintf('[Stage 7] Launching electrode naming GUI (%d clusters)...\n', N);

            % Load pial surface for display (use lh for visual context;
            % both hemispheres are plotted together)
            sumaDir = fullfile(self.locDirs.fs_subj, 'SUMA');
            lhFile  = fullfile(sumaDir, 'lh.pial.gii');
            rhFile  = fullfile(sumaDir, 'rh.pial.gii');
            lhSurf  = gifti(lhFile);
            rhSurf  = gifti(rhFile);

            % Open 3-D figure
            fig = figure('Name', sprintf('Electrode Naming — %s', self.subj), ...
                         'NumberTitle', 'off', 'Color', 'k');
            ax  = axes('Parent', fig, 'Color', 'k');
            hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');
            view(ax, 3);
            camlight(ax, 'headlight');
            material(ax, 'dull');

            % Plot pial surfaces
            patch(ax, 'Faces', lhSurf.faces, 'Vertices', lhSurf.vertices, ...
                'FaceColor', [0.75 0.70 0.65], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.4);
            patch(ax, 'Faces', rhSurf.faces, 'Vertices', rhSurf.vertices, ...
                'FaceColor', [0.75 0.70 0.65], 'EdgeColor', 'none', ...
                'FaceAlpha', 0.4);

            % Plot all clusters as grey dots with index labels
            scatter3(ax, self.clusters(:,1), self.clusters(:,2), self.clusters(:,3), ...
                60, repmat([0.5 0.5 0.5], N, 1), 'filled');
            for k = 1:N
                text(ax, self.clusters(k,1), self.clusters(k,2), self.clusters(k,3), ...
                    sprintf('  %d', k), 'Color', 'w', 'FontSize', 8);
            end
            title(ax, sprintf('%s — %d clusters detected', self.subj, N), ...
                'Color', 'w', 'FontSize', 11);

            % Initialise output table
            chanNames_out = cell(N, 1);
            types_out     = cell(N, 1);
            isArtifact    = false(N, 1);
            highlight     = [];   % handle to highlighted scatter point

            % Iterate over clusters
            k = 1;
            while k <= N

                % Highlight the current cluster
                delete(highlight);
                highlight = scatter3(ax, ...
                    self.clusters(k,1), self.clusters(k,2), self.clusters(k,3), ...
                    120, [1 0.8 0], 'filled', 'MarkerEdgeColor', 'w');
                drawnow;

                % Build dialog
                prompt    = {sprintf('Cluster %d of %d  (xyz: %.1f %.1f %.1f)\n\nChannel name:', ...
                                k, N, self.clusters(k,:))};
                dlgTitle  = 'Name this electrode';
                defaults  = {''};

                % If chanNames available, offer a listdlg first
                if ~isempty(self.chanNames)
                    remaining = setdiff(self.chanNames, chanNames_out, 'stable');
                    choices   = [remaining; {'[artifact — skip]'}; {'[free text]'}];
                    [sel, ok] = listdlg( ...
                        'ListString',   choices, ...
                        'SelectionMode','single', ...
                        'Name',         sprintf('Cluster %d / %d', k, N), ...
                        'PromptString', sprintf('xyz: %.1f  %.1f  %.1f', ...
                                                self.clusters(k,:)));
                    if ~ok
                        % User closed dialog — offer to go back or quit
                        choice = questdlg('Continue or go back?', '', ...
                            'Go back', 'Skip cluster', 'Quit naming', 'Go back');
                        switch choice
                            case 'Go back'
                                k = max(1, k-1);
                                continue;
                            case 'Skip cluster'
                                isArtifact(k) = true;
                                k = k + 1;
                                continue;
                            case 'Quit naming'
                                break;
                            otherwise
                                break;
                        end
                    end

                    selected = choices{sel};
                    if strcmp(selected, '[artifact — skip]')
                        isArtifact(k) = true;
                        k = k + 1;
                        continue;
                    elseif strcmp(selected, '[free text]')
                        answer = inputdlg(prompt, dlgTitle, 1, defaults);
                        if isempty(answer) || isempty(strtrim(answer{1}))
                            isArtifact(k) = true;
                            k = k + 1;
                            continue;
                        end
                        chanName = strtrim(answer{1});
                    else
                        chanName = selected;
                    end

                else
                    % No chanNames — free text only
                    choices = {'[artifact — skip]'; '[free text]'};
                    [sel, ok] = listdlg( ...
                        'ListString',   choices, ...
                        'SelectionMode','single', ...
                        'Name',         sprintf('Cluster %d / %d', k, N), ...
                        'PromptString', sprintf('xyz: %.1f  %.1f  %.1f', ...
                                                self.clusters(k,:)));
                    if ~ok || sel == 2
                        answer = inputdlg(prompt, dlgTitle, 1, defaults);
                        if isempty(answer) || isempty(strtrim(answer{1}))
                            isArtifact(k) = true;
                            k = k + 1;
                            continue;
                        end
                        chanName = strtrim(answer{1});
                    else
                        isArtifact(k) = true;
                        k = k + 1;
                        continue;
                    end
                end

                % Ask depth vs subdural
                typeChoice = questdlg( ...
                    sprintf('%s — electrode type?', chanName), ...
                    'Electrode type', 'depth', 'subdural', 'depth');
                if isempty(typeChoice)
                    typeChoice = 'depth';
                end

                chanNames_out{k} = chanName;
                types_out{k}     = typeChoice;

                % Update scatter colour: depth=blue, subdural=green
                if strcmp(typeChoice, 'depth')
                    col = [0.3 0.6 1.0];
                else
                    col = [0.3 1.0 0.5];
                end
                scatter3(ax, self.clusters(k,1), self.clusters(k,2), self.clusters(k,3), ...
                    60, col, 'filled');
                text(ax, self.clusters(k,1), self.clusters(k,2), self.clusters(k,3), ...
                    sprintf('  %s', chanName), 'Color', 'w', 'FontSize', 7);
                drawnow;

                k = k + 1;
            end

            delete(highlight);

            % Collect non-artifact entries
            valid     = ~isArtifact & ~cellfun(@isempty, chanNames_out);
            chanNames_out = chanNames_out(valid);
            types_out     = types_out(valid);
            xyz_valid     = self.clusters(valid, :);

            self.leads = table(chanNames_out, ...
                               xyz_valid(:,1), xyz_valid(:,2), xyz_valid(:,3), ...
                               types_out, ...
                               'VariableNames', {'chanName','x','y','z','type'});

            fprintf('[Stage 7] Naming complete: %d contacts accepted, %d marked as artifact.\n', ...
                sum(valid), sum(isArtifact));
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
