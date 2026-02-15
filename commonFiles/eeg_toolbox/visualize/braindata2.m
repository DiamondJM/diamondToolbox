% braindata2
%
% Responsibile for:
%   subject data (e.g. element_info)
%   tal (zLocalize) data
%   ROIs (subject + standard)
%   User data
%   Statistics
%   
% NOT responsible for plotting, geometries, etc
% Does NOT contain any references to BrainPlotter
%
%
% BrainPlotter responsible for:
%   Storing surface geometries
%   Data sets (underlays like curvature, sulc, atlases, + overlays like user data)
%   XYZ electrode points
%   
% braindata2 Properties:
%   rootEEGdir;     - subject's parent dir string
%   subj;           - subject string
%   filepaths;      - struct containing paths to files
%   userData;       - struct for all user data
%   roi;            - struct 
%   docs;           - struct for general subject data
%   tal;            - struct for tal data
%   zloc;           - struct for localization-algorithm data (e.g. CT-registered coords, anchors)
%   stats;          - struct for stats data
%   misc;           - miscellaneous struct
%
% braindata2 Methods:
%   braindata2      - Constructor. Pass subject and rootEEGdir to load
%   loadSubject     - Loads a subject
%   loadAverage     - Load fsaverage
%   set_default_filepaths - (sets filepaths property to defaults)
%   check_files     - (Display to user all loaded data, filepaths, and status)
%   load_all_files  - (Calls all the loading methods)
%   load_docs       - (Loads docs)
%   load_roi        - (Loads roi)
%   load_tal        - (Loads tal)
%   load_zloc       - (Loads zloc)
%   setROICs        - Sets number of ROICs if you need less than 2400
%   leadData2ROI    - Map channel data onto their ROI
%   leadData2ROIC   - Map channel data onto their ROIC
%   leadData2ROICMat- Map channel data onto their ROIC, matrix-based and faster
%   vertex2ROI      - Map a vertex to its ROICs/ROIs
%   roic2roi        - Map a ROIC to its ROI
%
% Static Methods:
%   xLoadSubj       - Load a list of braindata2 objects corresponding to a set of subjects
%   xLeadData2ROI   - Map channel data from a list of subjects onto their ROIC
%   peek            - Utility function; displays a preview of a table
%   
%
% For More help:
%       - Create an instance of brainplotter by calling "bd = braindata(subj, rootEEGdir)".
%       - You can use the MATLAB functions "methods(bd)" and "properties(bd)" to investigate your object.
%       - Try typing: "doc braindata2".
%       - See the tutorial on /Volumes/Shares/FRNU/labmeetings/2018_03_20_trotta_braindata/bdtutorial/braindata_tutorial.m
%         Please copy it locally before running it.  It has a FAQ and a walk-through. 
%
% See also: brainplotter, braindata

% REVISION HISTORY
%   ( ** Note: remember to update the version() function! ** )
%   04/2018 MST     - Use getToolboxDir instead of hard-coded path
%   05/2018 MST     - Change location of the glocal ROIC list
%   06/2018 MST     - getElementInfo instead of readtable
%                   - add leadData2ROICMat function
%   07/2018 MST     - chans2depth bug fixed

classdef braindata2 < handle
    
    % Constructor (and other standard class functions)
    methods 
        function bd = braindata2(subj, rootEEGdir, isSilent)
            % The input variables are completely optional, but
            % if you give them, it will save you a load_subject call
            if nargin < 3
                bd.isSilent = 0; 
            else
                bd.isSilent = isSilent;
            end
            
            if nargin >= 1, bd.subj = subj; end
            if nargin >= 2, bd.loadSubject(subj, rootEEGdir); end
            
           
            bd.set_default_filepaths();
            
            
        end
        
        function disp(bd)
            bd.line(0); bd.line(1); bd.line(0);
            fprintf('Subject    : %s\n', bd.subj);
            fprintf('RootEEGdir : %s\n', bd.rootEEGdir);
            fprintf('\n');
            
            disp('Tal:');
            disp(bd.tal);
            disp('ROI:');
            disp(bd.roi);
            disp('Docs:');
            disp(bd.docs);
            disp('Zloc:');
            disp(bd.zloc);
            fprintf('\n');
            bd.version();
            bd.line(0);
        end
    end
    
    
    % Public properties
    properties (Access = public) 
        
        rootEEGdir;     % subject's parent dir string
        subj;           % subject string
        filepaths;      % struct containing paths to files
        userData;       % struct for all user data
        roi;            % struct 
        docs;           % struct for general subject data
        tal;            % struct for tal data
        zloc;           % struct for localization-algorithm data (e.g. CT-registered coords, anchors)
        stats;          % struct for stats data
        misc;           % miscellaneous struct
        
        isSilent;       % whether braindata should be more silent in terms of its output
    end
    
    properties (Constant)
        stdNumVertices = 198812;
    end
    
    % Data loading, file reading methods
    methods (Access = public)
        
        function loadSubject(bd, subj, rootEEGdir)
            % Populates subject-specific docs, tal, roi, and zLoc data
            bd.subj = subj;
            bd.rootEEGdir = rootEEGdir;
            bd.set_default_filepaths();
            bd.check_files()
            bd.load_all_files();
        end
        function loadAverage(bd)
            bd.subj = 'fsaverage';
            bd.set_default_filepaths();
            
            bd.load_roi();
            
            if exist(bd.filepaths.roic_roi_lh, 'file')
                bd.roi.roic_roi_lh = readtableSafe(bd.filepaths.roic_roi_lh);
            end
            
            if exist(bd.filepaths.roic_roi_rh, 'file')
                bd.roi.roic_roi_rh = readtableSafe(bd.filepaths.roic_roi_rh);
            end
            
        end
        function set_default_filepaths(bd)
            % Does NOTHING except set default file paths
            
            bd.filepaths.visData = fullfile(getToolboxDir(), 'visualize/data');

            % Global ROI files
            bd.filepaths.ROIC_lh = fullfile(bd.filepaths.visData, 'zLocalize_ROI_centers_lh.mat');
            bd.filepaths.ROIC_rh = fullfile(bd.filepaths.visData, 'zLocalize_ROI_centers_rh.mat');
            
            if strcmpi(bd.subj, 'fsaverage')
                bd.rootEEGdir = bd.filepaths.visData;
                dsubj                          = fullfile(bd.rootEEGdir, bd.subj);
                bd.filepaths.surfaces          = fullfile(dsubj, 'SUMA');
                bd.filepaths.roic_roi_lh       = fullfile(dsubj, 'ROIC_mesh_LUT_lh.csv');
                bd.filepaths.roic_roi_rh       = fullfile(dsubj, 'ROIC_mesh_LUT_rh.csv');
                
            elseif ~isempty(bd.subj) && ~isempty(bd.rootEEGdir)
                dsubj   = fullfile(bd.rootEEGdir, bd.subj);

                bd.filepaths.element_info      = fullfile(dsubj, 'docs/element_info.csv');
                bd.filepaths.jacksheet         = fullfile(dsubj, 'docs/jacksheetMaster.csv');
                bd.filepaths.zloc              = fullfile(dsubj, 'tal/zloc');
                bd.filepaths.leads_xyz         = fullfile(dsubj, 'tal/leads.csv');
                bd.filepaths.leads_bp          = fullfile(dsubj, 'tal/zloc/mr_pre/coords_mid_euclid.csv');
                bd.filepaths.atlas             = fullfile(dsubj, 'tal/atlas/atlas_monopolar.csv');
                bd.filepaths.atlas_bp          = fullfile(dsubj, 'tal/atlas/atlas_bipolar.csv');
                bd.filepaths.atlas_sim         = fullfile(dsubj, 'tal/atlas/atlas_simple.csv');
                bd.filepaths.atlas_sim_bp      = fullfile(dsubj, 'tal/atlas/atlas_bipolar_simple.csv');
                bd.filepaths.ct_xyz            = fullfile(dsubj, 'tal/zloc/CT_1/transform/coords.csv');
                bd.filepaths.ct_2_xyz          = fullfile(dsubj, 'tal/zloc/CT_2/transform/coords.csv');
                bd.filepaths.roic_lead_mono    = fullfile(dsubj, 'tal/roi/lead_ROIC_LUT_monopolar.mat');
                bd.filepaths.roic_lead_bipo    = fullfile(dsubj, 'tal/roi/lead_ROIC_LUT_bipolar.mat');
                
                % bd.filepaths.roic_roi_lh       = fullfile('/Users/diamondjm/Documents/Josh_work/commonFiles/Braindata',bd.subj,'ROIC_mesh_LUT_lh.csv');
                % bd.filepaths.roic_roi_rh       = fullfile('/Users/diamondjm/Documents/Josh_work/commonFiles/Braindata',bd.subj,'ROIC_mesh_LUT_rh.csv');
                bd.filepaths.roic_roi_lh       = fullfile(dsubj, 'tal/roi/ROIC_mesh_LUT_lh.csv');
                bd.filepaths.roic_roi_rh       = fullfile(dsubj, 'tal/roi/ROIC_mesh_LUT_rh.csv');
                
                bd.filepaths.lead_mesh         = fullfile(dsubj, 'tal/roi/lead_mesh_LUT_monopolar.mat');
                bd.filepaths.lead_mesh_bp      = fullfile(dsubj, 'tal/roi/lead_mesh_LUT_bipolar.mat');
                bd.filepaths.surfaces          = fullfile(dsubj, 'tal/zloc/freesurfer', bd.subj, 'SUMA');
                bd.filepaths.anchors           = fullfile(dsubj, 'tal/zloc/anchors/anchors.csv');
                bd.filepaths.snap_1_xyz        = fullfile(dsubj, 'tal/zloc/mr_pre/coords_snap_1.csv');
                bd.filepaths.snap_2_xyz        = fullfile(dsubj, 'tal/zloc/mr_pre/coords_snap_2.csv');
                bd.filepaths.proj_1_xyz        = fullfile(dsubj, 'tal/zloc/mr_pre/coords_1.csv'); 
                bd.filepaths.proj_2_xyz        = fullfile(dsubj, 'tal/zloc/mr_pre/coords_2.csv'); 
                bd.filepaths.dural_dist_1      = fullfile(dsubj, 'tal/zloc/mr_pre/dural_dist_1.csv'); 
                bd.filepaths.dural_dist_2      = fullfile(dsubj, 'tal/zloc/mr_pre/dural_dist_2.csv'); 
                bd.filepaths.utah_xyz          = fullfile(dsubj, 'tal/zloc/anchors/utah.csv'); 
                
            end
            
            
            % Shouldn't trigger alarm if these are missing
            bd.misc.optional_filepaths = {...
                'snap_2_xyz'
                'dural_dist_2'
                'proj_2_xyz'
                'ct_2_xyz'
                'utah_xyz'};
        end

        function check_files(bd)
            % Displays file status
            if bd.isSilent
                return
            end
            bd.printf('Checking files....\n');
            bd.line(1);
            fn = fieldnames(bd.filepaths);
            bd.printf('\n%15s\t%12s\t%25s\t%s\n\n','Property Name', '(!) Status', 'Filename', 'Full Path');
            for i = 1:numel(fn)
                x = exist(bd.filepaths.(fn{i}), 'file');
                if x > 0
                    status = 'EXISTS';
                elseif ismember(fn{i}, bd.misc.optional_filepaths)
                    status = '   -   ';
                else
                    status = '! MISSING';
                end
                [~,fname,ext] = fileparts(bd.filepaths.(fn{i}));
                bd.printf('%15s\t%12s\t%25s\t%s\n', fn{i}, status, [fname ext], bd.filepaths.(fn{i}));
            end
            bd.line(0);
        end
        function load_all_files(bd)
            % Loads all filepaths
            bd.load_docs();
            bd.load_roi();
            bd.load_tal();
            bd.load_zloc(); 
            if ~bd.isSilent
                bd.disp();
            end

            
        end % load_all_files
        function load_docs(bd)
            bd.printf('Loading docs....'); 
            if exist(bd.filepaths.jacksheet, 'file')
                bd.docs.jacksheet = readtableSafe(bd.filepaths.jacksheet);
            end
            
            % if exist(bd.filepaths.element_info, 'file')
            %     bd.docs.element_info = getElementInfo([],[], bd.filepaths.element_info);
            % end
            bd.printf('Done!\n');
        end
        function load_roi(bd)
            bd.printf('Loading roi....');
            bd.roi.n = 2400;
            
            if isfield(bd.filepaths, 'roic_lead_mono') && exist(bd.filepaths.roic_lead_mono, 'file')
                data = load(bd.filepaths.roic_lead_mono);
                bd.roi.roic_lead_mono = data.(char(fieldnames(data)));
            end
            
            if isfield(bd.filepaths, 'roic_lead_bipo') && exist(bd.filepaths.roic_lead_bipo, 'file')
                data = load(bd.filepaths.roic_lead_bipo);
                bd.roi.roic_lead_bipo = data.(char(fieldnames(data)));
            end

            localizer_rois(bd.subj,bd.rootEEGdir,1);
            % JD edit: force new localizer to run. 
            fprintf('There was an issue with how the localization data was originally run, where only a subset of ROIs were retained. Therefore we have to re-run the ROI generation procedure at this time.'); 
            % In the future, the whole localization pipeline could be rerun
            % with the proper change made to localizer_rois. In that case,
            % this call would not be necessary. 

            
            if isfield(bd.filepaths, 'roic_roi_lh') && exist(bd.filepaths.roic_roi_lh, 'file')
                bd.roi.roic_roi_lh = readtableSafe(bd.filepaths.roic_roi_lh);
            end
            
            if isfield(bd.filepaths, 'roic_roi_rh') && exist(bd.filepaths.roic_roi_rh, 'file')
                bd.roi.roic_roi_rh = readtableSafe(bd.filepaths.roic_roi_rh);
            end
            
            if isfield(bd.filepaths, 'ROIC_lh') && exist(bd.filepaths.ROIC_lh, 'file')
                data = load(bd.filepaths.ROIC_lh);
                bd.roi.ROIC_lh = data.vert_idx;
            end
            
            if isfield(bd.filepaths, 'ROIC_rh') && exist(bd.filepaths.ROIC_rh, 'file')
                data = load(bd.filepaths.ROIC_rh);
                bd.roi.ROIC_rh = data.vert_idx;
            end
            
            if isfield(bd.filepaths, 'lead_mesh') && exist(bd.filepaths.lead_mesh, 'file')
                data = load(bd.filepaths.lead_mesh);
                bd.roi.lead_mesh = data.lead_mesh_lut;
            end
            
            if isfield(bd.filepaths, 'lead_mesh_bp') && exist(bd.filepaths.lead_mesh_bp, 'file')
                data = load(bd.filepaths.lead_mesh_bp);
                bd.roi.lead_mesh_bp = data.lead_mesh_lut;
            end
            
            bd.printf('Done!\n');
        end
        function load_tal(bd)
            bd.printf('Loading tal....'); 
            if exist(bd.filepaths.leads_xyz, 'file')
                bd.tal.xyz = readtableSafe(bd.filepaths.leads_xyz);
            end
            
            if exist(bd.filepaths.leads_bp, 'file')
                bd.tal.xyz_bp = readtableSafe(bd.filepaths.leads_bp);
            end
            
            if exist(bd.filepaths.atlas, 'file')
                bd.tal.atlas = readtableSafe(bd.filepaths.atlas);
            end
            
            if exist(bd.filepaths.atlas_bp, 'file')
                bd.tal.atlas_bp = readtableSafe(bd.filepaths.atlas_bp);
            end
            
            if exist(bd.filepaths.atlas_sim, 'file')
                bd.tal.atlas_sim = readtableSafe(bd.filepaths.atlas_sim);
            end
            
            if exist(bd.filepaths.atlas_sim_bp, 'file')
                bd.tal.atlas_sim_bp = readtableSafe(bd.filepaths.atlas_sim_bp);
            end
            
            if exist(bd.filepaths.utah_xyz, 'file')
                bd.tal.xyz_utah = readtableSafe(bd.filepaths.utah_xyz);
            end
            
            
            bd.printf('Done!\n');
        end
        function load_zloc(bd)
            bd.printf('Loading zloc....'); 
            if exist(bd.filepaths.ct_xyz, 'file')
                bd.zloc.ct_xyz = readtableSafe(bd.filepaths.ct_xyz);
            end
            
            if exist(bd.filepaths.ct_2_xyz, 'file')
                bd.zloc.ct_2_xyz = readtableSafe(bd.filepaths.ct_2_xyz);
            end
            
            if exist(bd.filepaths.anchors, 'file')
                bd.zloc.anchors = readtableSafe(bd.filepaths.anchors);
            end
            
            if exist(bd.filepaths.snap_1_xyz, 'file')
                bd.zloc.snap_1_xyz = readtableSafe(bd.filepaths.snap_1_xyz);
            end
            
            if exist(bd.filepaths.snap_2_xyz, 'file')
                bd.zloc.snap_2_xyz = readtableSafe(bd.filepaths.snap_2_xyz);
            end
            
            if exist(bd.filepaths.dural_dist_1, 'file')
                bd.zloc.dural_dist_1 = readtableSafe(bd.filepaths.dural_dist_1);
            end
            
            if exist(bd.filepaths.dural_dist_2, 'file')
                bd.zloc.dural_dist_2 = readtableSafe(bd.filepaths.dural_dist_2);
            end
            
            if exist(bd.filepaths.proj_1_xyz, 'file')
                bd.zloc.proj_1_xyz = readtableSafe(bd.filepaths.proj_1_xyz);
            end
            
            if exist(bd.filepaths.proj_2_xyz, 'file')
                bd.zloc.proj_2_xyz = readtableSafe(bd.filepaths.proj_2_xyz);
            end
            bd.printf('Done!\n');
        end
        
    end
    
    % ROI
    methods 
        function setROICs(bd, n)
            % setROICs(n) sets the number of ROIs to use. (max is 2400)
            % bd.roic_roi table is filtered down to the subset 1:n
            if nargin < 2, n=2400; end
            
            % User-determined for how many ROIs
            bd.roi.n = n;
            
            % filter table to first n ROIs
            try
                bd.roi.roic_roi_lh = bd.roi.roic_roi_lh(bd.roi.roic_roi_lh.ROIC_ndx <= n, :);
                bd.roi.roic_roi_rh = bd.roi.roic_roi_rh(bd.roi.roic_roi_rh.ROIC_ndx <= n, :);
            catch
            end

        end
        function [VPerROI, valPerROI, roiMask, rh_begin, VFind] = leadData2ROI(bd, data, chanNames, varargin)
            % [VPerROI, valPerROIC, sparseMask] = leadData2ROI(bd, data, chanNames, ...)
            %
            %   Averages channel data at each ROI center
            %
            % INPUT
            %   data        - 1 x n array of numeric data per channel
            %   chanNames   - 1 x n cell array of channels to plot (use a-b for bipolars)
            %
            % OPTIONAL KEY-VALUE INPUT
            %   rroi        - size of ROI (mm). Default 9
            %   force_bilat - If 1, return [lh; rh] regardless of channels given (good for cross-subj). Default 0.
            %   agg_fun     - Function to apply to channel data at each ROIC. Default @nanmean
            %   VPerROIMat  - Pass in adjacency matrix from VPerROIMat if you have it.
            %   Mat         - If 1, mat version. VPerROI will be adjecency matrix instead of cells.
            %   VFind
            %
            % OUTPUT
            %   VPerROI     - Mesh vertices at all ROICs with data
            %   valPerROIC  - mean data value at all ROICs with data aggregated from given channel data
            %   sparseMask  - logical index into bd.roi.ROIC_[l/r]h.
            %   rh_begin    - index at which right hemisphere values/vertices begin
            %   VFind       - 
            %
            % NOTES
            %   * sum(sparseMask) == size(VperROI,1) == size(valPerROI,1) == {one of bd.roi.n or 2*bd.roi.n}
            %   * If you give channels from both hemispheres, you will get twice as many resuts (2 * bd.roi.n), returned
            %     as [left; right]. This is also the case if you force_bilat
            %   * In the bilateral case, you will need rh_begin to split your VPerROI and valPerROI results
            %
            % See Also: braindata2/leadData2ROIC
            
            % input
            ip = inputParser; 
            ip.KeepUnmatched = 1;
            ip.addParameter('rroi',9);
            ip.addParameter('force_bilat',0);
            ip.addParameter('agg_fun', @nanmean); % @nanmean
            ip.addParameter('VPerROIMat', []);
            ip.addParameter('Mat', 0);
            ip.addParameter('VFind', []);
            ip.parse(varargin{:});
            rroi = ip.Results.rroi;
            force_bilat = ip.Results.force_bilat;
            agg_fun = ip.Results.agg_fun;
            VPerROIMat = ip.Results.VPerROIMat;
            use_mat = ip.Results.Mat;
            VFind = ip.Results.VFind;
            
            % output
            VPerROI     = [];
            valPerROI   = [];
            roiMask     = [];
            rh_begin    = 0;

            
            assert(isfield(bd.docs,'jacksheet'), 'docs.jacksheet must be loaded first')
            
            % Map lead data to ROIC
            if use_mat
                valPerROI = full(bd.leadData2ROICMat(data, chanNames, varargin{:}));
                roiMask(find(valPerROI)) = true;
                
                
            else
                [valsPerROIC, roiMask] = bd.leadData2ROIC(data, chanNames, varargin{:});
                
                % Aggregate data at ROIC
                ndxs = find(roiMask);
                for i=1:sum(roiMask)
                    valPerROI(i) = agg_fun(valsPerROIC(ndxs(i),:));
                end
            end
            
            if ~use_mat || isempty(VPerROIMat) || isempty(VFind)
                % calculate VPerROIMat
                % Find ROI vertex set per ROIC
                [left, right] = bd.chans2hem(chanNames, bd.docs.jacksheet);
                if ~isempty(left) %|| force_bilat
                    VPerROI = bd.roic2roi(bd.roi.ROIC_lh(1==roiMask(1:bd.roi.n)), 'lh', rroi);
                end
                if ~isempty(right) %|| force_bilat
                    VPerROI = [VPerROI; bd.roic2roi(bd.roi.ROIC_rh(roiMask(end-bd.roi.n+1:end)), 'rh', rroi)];
                    rh_begin = 1 + length(find(roiMask(1:bd.roi.n)));
                elseif force_bilat
                    rh_begin = 1 + length(find(roiMask(1:bd.roi.n)));
                else
                    rh_begin = 0;
                end

                if isempty(left), rh_begin = 1; end

                if use_mat
                    [VPerROI,~,VFind] = bd.VPerROI2Mat(VPerROI, roiMask);
                    roiMask = logical(roiMask);
                    VPerROI = VPerROI(roiMask);
                end
                
            else
                VPerROI = VPerROIMat;
            end
            
        end
        
        function [VPerROIMat, vertexCnt, VFind] = VPerROI2Mat(bd, VPerROI, roiMask, useSparse)
            if nargin < 4, useSparse = 0; end
            
            % VPerROIMat = VPerROI2Mat(bp, VPerROI) is a helper function used by plotRegionsDataMat to
            % convert a cell array to a Vertex X ROIC adjacency matrix, normalized
            VPerROIMat = zeros(bd.stdNumVertices, numel(roiMask));
            roiFind = find(roiMask);
            for i = 1:numel(roiFind)
                VPerROIMat(VPerROI{i}, roiFind(i)) = 1; 
            end
            
            % We want matrix multiplication (dot product) to, at each vertex, average all ROI values to which the 
            % vertex belongs. Therefore, normalize the VPerROIMat rows.
            if useSparse
                VPerROIMat = sparse(VPerROIMat);
            end
            
            vertexCnt = sum(VPerROIMat, 2);
            VFind = find(vertexCnt);
            VPerROIMat(VFind,:) = VPerROIMat(VFind,:) ./ vertexCnt(VFind);
            VPerROIMat = VPerROIMat(VFind,:);
            
            
        end
        
        function [valsPerROICMat, chanByROICMat] = leadData2ROICMat(bd, dataMat, chanNames, varargin)
            % loadData2ROICmat maps an N X Electrode data matrix to an N X ROIC           
            %   
            % INPUT
            %   dataMat     - n x c array of numeric data per channel.
            %   chanNames   - cell array of length c of channels corresponding to dataMat.
            %
            % OPTIONAL KEY-VALUE INPUT
            %   rroi        - size of ROI (mm). Default 9
            %   force_bilat - If 1, return [lh, rh] ROIs regardless of channels given (good for cross-subj). Default 0.
            %   useMean     - If 1 (default), take mean of electrode data per ROI. If 0, take sum.
            %
            % OUTPUT
            %   valsPerROICMat - sparse matrix of N x ROI
            %   chanByROICMat  - this corresponds to GIVEN chanNames (note this is different from the bd property)
            %
            %
            % If you pass force_bilar or subject has a bilateral implant, you will get twice as many
            % rows, where rows 1:n are left hemisphere, and rows n+1:2n are right hemisphere ROIs, where
            % n is bd.roi.n
            %
            % NOTES:
            %   This function is faster than leadData2ROIC, relying on (sparse) matrix
            %   multiplication:
            %
            %   [ n X electrode ] * [ electrode X ROIC ]
            %       ^ dataMat           ^ ROICs within D mm
            %
            %   The result is an N X ROIC matrix where each element represents
            %   either:
            %       1) The electrode-average at each ROI   OR
            %       2) The elecrode-cumulative at each ROI (if ~useMean)
            %
            %   Note that the dataMat has arbitrary units, but the electrode and
            %   ROI matrices are binary (similar to adjacency matrices).
            %   Prior to multiplying the dataMat,
            %   if we are going to *average* at each ROI (instead of cumulative)
            %   we make sure each column of the [electrode X ROIC] matrix is scaled
            %   so that its elements sum to one, thereby taking the average.
            %   
            %
            % See Also: braindata2/leadData2ROIC

            % parse input
            ip = inputParser;
            ip.KeepUnmatched = 1;
            ip.addParameter('rroi',9);
            ip.addParameter('force_bilat',0);
            ip.addParameter('useMean', 1);
            ip.parse(varargin{:});
            rroi = ip.Results.rroi;
            force_bilat = ip.Results.force_bilat;
            useMean = ip.Results.useMean;

            if size(dataMat,1) > 1 && size(dataMat,2) == 1
                dataMat = dataMat';
            end
            assert(size(dataMat, 2) == numel(chanNames), 'chanNames length (%d) must match size(dataMat, 2) (%d)', numel(chanNames), size(dataMat,2));
            %isBipolar = all(cellfun(@(x) strfound(x,'-'), chanNames));
            %assert(~isBipolar, 'leadData2ROICMat not implemented for bipolars');
            
            % Build/store per hemisphere (store result in object)
            chans_all = bd.roi.roic_lead_mono.chanName;

            % At this point, we have the mat's built for all channels ordered in the same order as
            % bd.roi.roic_lead_mono. We now have to order by the given channel names and do the multiplication
            missing = setdiff(chanNames, chans_all);
            if ~isempty(missing)
                % This is a serious error that indicates incorrect roi data files
                disp('\nThe following channels are missing from SUBDURALS:');
                disp(missing);
                error('bd.roi.roic_lead_mono is missing the above channels. Either you gave incorrect channels OR localizer_rois.m needs to be re-run!');
            end
            
            % Retrieve or Calculate/store the chanByROICmat for each hemisphere
            if isfield(bd.roi, 'chanByROICmat')
                chanByROICMat = bd.roi.chanByROICmat;
            else
                chanByROICMat = bd.buildMat(rroi);
                bd.roi.chanByROICmat = chanByROICMat;
            end
            
            % These arrays correspond to the *given* chanNames and
            % provide us the indices into the saved chans (and save data)
            [~, saved2given] = ismember(chanNames, chans_all);
            chanByROICMat = chanByROICMat(saved2given, :);
            
            % alter mat so that....
            % columns entries should sum to one so that the dot
            % product of a column with a data vector results in an
            % average of electrodes in an ROI
            if useMean
                n_elec_per_roi = full(sum(chanByROICMat,1));
                mask = find(n_elec_per_roi);
                chanByROICMat(:,mask) = chanByROICMat(:,mask) ./ n_elec_per_roi(mask);
            end

            % Transform
            valsPerROICMat = dataMat * chanByROICMat;
            
            if ~force_bilat
                [chans_lh, chans_rh] = bd.chans2hem(chans_all, bd.docs.jacksheet);
                use_lh = ~isempty(intersect(chans_lh, chanNames));
                use_rh = ~isempty(intersect(chans_rh, chanNames));
                if ~use_lh
                    valsPerROICMat = valsPerROICMat(:,bd.roi.n + 1 : end);
                elseif ~use_rh
                    valsPerROICMat = valsPerROICMat(:,1 : bd.roi.n);
                end
            end
        end
        
        function mat = buildMat(bd, rroi)
            % Private function called by leadData2ROICMat
            % ---------------------------------------------------
            % Construct electrode -> ROI map (elec X ROIC matrix)
            % ---------------------------------------------------
            %   This is done by building a sparse matrix from the data read
            %   in the existing electrode -> ROIC_mesh_ndx table, [lh rh]
            % note roics are the mesh vertices in a hemisphere.
            % note lh comes first
            i = [];
            j = [];
            t = bd.roi.roic_lead_mono;
            m = height(t);
            for e = 1:m
                % for each electrode, make a row in the matrix...
                V = t.ROIC_mesh_ndx{e};        % ROIC vertices within D
                V = V(t.ROIC_d{e} <= rroi);

                if strcmpi(t.whichHemi{e}, 'lh');
                    roics = bd.roi.ROIC_lh;
                    eroic = arrayfun(@(v) find(roics == v, 1), V); % ROIC index in ROIC list
                    
                elseif strcmpi(t.whichHemi{e}, 'rh');
                    roics = bd.roi.ROIC_rh;
                    eroic = bd.roi.n + arrayfun(@(v) find(roics == v, 1), V); % ROIC index in ROIC list + 2400
                    
                else
                    error('Electrode %s has no hemisphere marked!', t.chanName{e});
                end

                % for each ROI to which this electrode belongs, mark an entry for the sparse matrix
                cnt = numel(eroic);
                i = cat(1, i, repmat(e, cnt, 1));   % row index, add electrode index-into-chanNames
                j = cat(1, j, eroic);               % column index (add roic index)
                
            end
            mat = sparse(i, j, 1, m, 2*bd.roi.n);
            
        end
        function [valsPerROIC, sparseMask] = leadData2ROIC(bd, data, chanNames, varargin)
            % [valsPerROIC, sparseMask] = leadData2ROIC(bd, data, chanNames, ...)
            %
            %   Maps channel data onto their ROICs.
            %   This method is called by leadData2ROI, so if you are not sure which to use, you probably want
            %     to use leadData2ROI, *not this one*
            %
            % INPUT
            %   data        - 1 x n array of numeric data per channel
            %   chanNames   - 1 x n cell array of channels to plot (use a-b for bipolars)
            %
            % OPTIONAL KEY-VALUE INPUT
            %   rroi        - size of ROI (mm). Default 9
            %   force_bilat - If 1, return [lh; rh] regardless of channels given (good for cross-subj). Default 0.
            % 
            % OUTPUT
            %   valsPerROIC - matrix of ROIC x channel
            %   sparseMask  - mask into valsPerROIC rows which have any data
            %
            % data should be a vector with 1 value per bd.leads or you can pass specific chanNames
            % valsPerROIC has columns ordered by chanNames, rows ordered by ROIC
            % to get count in region, use sum(~isnan(roi_data),2);
            %
            % If you passed channel names that come from BOTH hemispheres, you will get twice as many
            % rows, where rows 1:n are left hemisphere, and rows n+1:2n are right hemisphere ROIs, where
            % n is bd.roi.n
            %
            % example usage:
            %   [valsPerROIC, sparseMask] = bd.data2ROIC(cgamma, chans, 0, RROI);
            %   roicd = nanmean(rd, 2);
            %   V = bd.roic2roi(bd.roi.roics(sparseMask), RROI);
            %   bd.plotRegionsData(V, roicd(mask), [], [], 1.0, 0);
            %
            % See also: braindata2/leadData2ROI, braindata2/leadData2ROICMat
                        
            ip = inputParser;
            ip.KeepUnmatched = 1;
            ip.addParameter('rroi',9);
            ip.addParameter('force_bilat',0);
            ip.parse(varargin{:});
            rroi = ip.Results.rroi;
            force_bilat = ip.Results.force_bilat;
            
            valsPerROIC = nan(bd.roi.n, 1);
            sparseMask  = false(bd.roi.n, 1);
            
            if numel(chanNames) == 0, return; end
            
            assert(numel(data) >= numel(chanNames));
            isBipolar = all(cellfun(@(x) strfound(x,'-'), chanNames));
            [left, right, li, ri] = bd.chans2hem(chanNames, bd.docs.jacksheet);
            
            % Handle each hemisphere separately and combine the results
            if ~isempty(left) && ~isempty(right) || force_bilat
                modargin = [varargin, {'force_bilat',0}]; 
                [valsPerROIC_l, sparseMask_l] = leadData2ROIC(bd, data(li), left, modargin{:});
                [valsPerROIC_r, sparseMask_r] = leadData2ROIC(bd, data(ri), right, modargin{:});
                wl = size(valsPerROIC_l, 2);
                wr = size(valsPerROIC_r, 2);
                w = max(wl, wr);
                h = size(valsPerROIC_l, 1);
                valsPerROIC_l = [valsPerROIC_l, nan(h, w - wl)];
                valsPerROIC_r = [valsPerROIC_r, nan(h, w - wr)];
                valsPerROIC = [valsPerROIC_l; valsPerROIC_r];
                sparseMask = logical([sparseMask_l; sparseMask_r]);
                return
                
            % These are the single-hemisphere case
            elseif ~isempty(left)
                hemisphere = 'lh';
                roic_hem = bd.roi.ROIC_lh;
            elseif ~isempty(right)
                hemisphere = 'rh';
                roic_hem = bd.roi.ROIC_rh;
            end
            
            if isBipolar
                leads = bd.roi.roic_lead_bipo.chanName;
            else
                leads = bd.roi.roic_lead_mono.chanName;
                assert(isempty(setdiff(chanNames, leads)), 'you gave chanNames not in leads');
            end

            % given chanNames index into leads; thus i'th data point must correspond to i'th chanNdx
            chanNdx = cellfun(@(chan) find(strcmpi(chan,leads),1), chanNames);   
            
            % for every roic we may have multiple data points
            valsPerROIC = nan(bd.roi.n, numel(data));
            
            if isBipolar
                std_nodes = bd.roi.roic_lead_bipo.nearest_mesh_ind;
            else
                std_nodes = bd.roi.roic_lead_mono.nearest_mesh_ind;
            end
              
            for i = 1:length(chanNames)
                % 1) get nearest pial mesh notes
                chan_nodes = std_nodes{chanNdx(i)};
                
                
                [roic_chan, ~] = bd.vertex2ROI(chan_nodes, hemisphere, rroi);
                roicNdx = find(ismember(roic_hem, roic_chan));
                valsPerROIC(roicNdx,i) = data(i);
            end
            % The above loop results in a sparse matrix where data is put into columns of each ROIC row
            
            
            % mask of any ROICs with datacount > 0;
            sparseMask = sum(~isnan(valsPerROIC),2) > 0;
        end
        function [roics, roi, joshOutput] = vertex2ROI(bd, vertex, hemisphere, rroi)
            % [roics, roi, roicNdx] = vertex2ROI(bd, vertex, rroi)
            % Maps a vertex to its ROICs and their ROIs
            % 
            %
            % Input:
            %   vertex      - vertex set to map
            %   hemisphere  - 'lh' or 'rh' (may omit if only 1 hemisphere)
            %
            % Input (Optional):
            %   rroi     - ROI radius (default Infinite)
            %
            % Output:
            %   roics       - roi centers as surface vertices
            %   roi         - vertices within rroi_mm of *any* ROI center (you might not want this*)
            %   roicNdx     - index into bd.roi.roics of roi centers
            
            
            % general procedure is vertices -----> ROIC -----> ROI vertices
            %                                 (1)         (2)
            
            % *Note: you might not want the roi output, because it isn't segregated into ROICs; it includes
            %       the union of all ROIs which contained any of the given vertices
                  
            if nargin < 3 || isempty(hemisphere)
                if      (~isfield(bd.roi,'ROIC_lh') || isempty(bd.roi.ROIC_lh)) && (isfield(bd.roi,'ROIC_rh') && ~isempty(bd.roi.ROIC_rh))
                    hemisphere = 'rh';
                elseif  (~isfield(bd.roi,'ROIC_rh') || isempty(bd.roi.ROIC_rh)) && (isfield(bd.roi,'ROIC_lh') && ~isempty(bd.roi.ROIC_lh))
                    hemisphere = 'lh';
                else
                    error('You must pass a hemisphere parameter');
                end
            end
            if nargin < 4 || isempty(rroi), rroi = Inf; end
            

            roic_roi_hem = sprintf('roic_roi_%s', hemisphere);
            roic_hem     = sprintf('ROIC_%s', hemisphere);
            assert(isfield(bd.roi, roic_roi_hem) && ~isempty(bd.roi.(roic_roi_hem)), 'Empty or non-existant, try reloading: %s', roic_roi_hem);
            
            t = bd.roi.(roic_roi_hem);   % note v, d, r columns: vertex, distance, center
            d_mask  = t.d <= rroi;
            v_mask  = ismember(t.vertex, vertex);  
            roics   = unique(t.ROIC_mesh_ndx(d_mask & v_mask),'stable');       % 1 roics which are within distance of query vertices
            
            %%%%%%%
            % JD edits: 
            %%%%%%%
            
            %roicNdx = unique(t.ROIC_ndx(d_mask & v_mask), 'stable');
            %assert(all(vector(roics) == vector(bd.roi.(roic_hem)(roicNdx))))
            c_mask = ismember(t.ROIC_mesh_ndx, roics);
            roi = unique(t.vertex(d_mask & c_mask));                            % 2 vertices which are within distance of centers
            %roicNdx = arrayfun(@(x) find(x == bd.roi.roics), roics,'uniformOutput',0);
            %roicNdx = cell2mat(roicNdx(~cellfun('isempty',roicNdx)));
            
            joshOutput = t(c_mask,:);
            
        end
        function VPerROI = roic2roi(bd, roics, hemisphere, rroi_mm)
            % returns sets of mesh indices which fall within rroi_mm of given roics
            %
            % INPUT
            %   roics       - n-length list of ROICs as *index-in-mesh*
            %   hemisphere  - lh or rh (may be omitted in unilateral implants
            %   rroi_mm     - (optional) geodesic radius to grow ROIs (usually subject-specific)
            %   
            %
            % OUTPUT
            %   vertexSets - n-length cell array of mesh indices for each given ROIC

            if nargin < 3 || isempty(hemisphere)
                if      (~isfield(bd.roi,'ROIC_lh') || isempty(bd.roi.ROIC_lh)) && (isfield(bd.roi,'ROIC_rh') && ~isempty(bd.roi.ROIC_rh))
                    hemisphere = 'rh';
                elseif  (~isfield(bd.roi,'ROIC_rh') || isempty(bd.roi.ROIC_rh)) && (isfield(bd.roi,'ROIC_lh') && ~isempty(bd.roi.ROIC_lh))
                    hemisphere = 'lh';
                else
                    error('You must pass a hemisphere parameter');
                end
            end
            if nargin < 4, rroi_mm = 9; end
            
            fieldname = sprintf('roic_roi_%s', hemisphere);
            assert(isfield(bd.roi, fieldname) && ~isempty(bd.roi.(fieldname)), 'Empty or non-existant roi field, try reloading: %s', fieldname);
            
            VPerROI = cell(length(roics), 1);

            % Use roi LUT to find ROI vertex set
            t = bd.roi.(fieldname);
            d_mask = (t.d <= rroi_mm);
            t = t(d_mask, :);
            for i = 1:length(roics)
                roic_mask = (t.ROIC_mesh_ndx == roics(i));
                roi_verts = unique(t{roic_mask, 'vertex'});
                VPerROI{i} = roi_verts;
            end
        end
    end
    
    % Private
    methods (Access='private')
        function printf(bd, varargin)
            if ~bd.isSilent
                fprintf(varargin{:});
            end
        end
        
        function line(bd, title)
            if nargin < 1, title=0; end
            if title
                bd.printf('================================================ Brain Data 2 ==================================================\n')
            else
                bd.printf('================================================================================================================\n')
            end
        end
    end
    
    % EZ "easy" methods
    % These methods are designed to reproduce braindata 1.0's functionality and demostrate how
    % a user should use braindata2 in conjunction with brainplotter
    methods
        
        function varargout = ezplot(bd, bp, hax)
            % ezplot is an easy way to make a figure and plot the current subject's brain
            % Look at the code to see how braindata2 and brainplotter are used
            %
            % Input:
            %   bp      - brainplotter 
            %   hax     - axis handle
            % 
            % Output:
            %   first optional - bp brainplotter object
            
            if nargin < 2 || isempty(bp), bp = []; end
            if nargin < 3 || isempty(hax), hax = []; end
            
            
            
            if isempty(hax)
                f = figure();
                axis();
                hax = f.CurrentAxes;
            else
                axes(hax); %- hAx passed might not be the current axis... so do it here. 
            end
            if isempty(bp) || isempty(fieldnames(bp.surfaces))
                % Create a brain plotter and load relevant hemisphere's surfaces
                bp = bd.ez_get_plotter();
            end
            
            % Plot surface
            bp.clear();
            
            surfs = fieldnames(bp.surfaces);
            for i=1:numel(surfs)
                bp.plot(surfs{i});
            end
            
            if nargout >= 1, varargout{1} = bp; end    
        end
        function varargout = ezplotEcog(bd, bp, hax)
            % ezplotEcog is an easy way to make a figure and plot the current subject's brain with Ecog electrodes
            % Look at the code to see how braindata2 and brainplotter are used
            %
            % Input:
            %   bp      - brainplotter 
            %   hax     - axis handle
            % 
            % Output:
            %   first optional - bp brainplotter object
            
            if nargin < 2 || isempty(bp), bp = []; end
            if nargin < 3 || isempty(hax), hax = []; end
            
            
            
            if isempty(hax)
                f = figure();
                axis();
                hax = f.CurrentAxes;
            end
            if isempty(bp) || isempty(fieldnames(bp.surfaces))
                % Create a brain plotter and load relevant hemisphere's surfaces
                bp = bd.ez_get_plotter();
            end
            
            % Plot surface
            bp.unplot();
            bp.clear();
            
            surfs = fieldnames(bp.surfaces);
            for i=1:numel(surfs)
                bp.plot(surfs{i});
            end
            
            % Plot all the subdural electrodes
            jack = bd.docs.jacksheet;
            subdurals = jack.chanName(strcmpi(jack.hardwareType, 'SUBDURAL'));
            xyz = bd.zloc.snap_1_xyz(ismember(bd.zloc.snap_1_xyz.chanName, subdurals),:);
            bp.plotPointGroup(xyz, 'element_info',bd.docs.element_info);
            
            if nargout >= 1, varargout{1} = bp; end   
            
            bp.legend('orientation', 'horizontal');
            set(gca,'FontSize', 16);
            
        end
        function bp = ez_get_plotter(bd)
            % private function for creating a brainplotter object and loading a pial surface
            bp = brainplotter();
                
            % Load pial version of each hemisphere in jacksheet
            if isfield('docs',bd) && isfield(bd.docs, 'jacksheet')
                hemis = unique(bd.docs.jacksheet.whichHemi);
                hemis = hemis(~cellfun(@isempty, hemis));
                for i = 1:numel(hemis)
                    bp.loadSurface(bd, ['pial ' hemis{i}]);
                end
            else
                bp.loadSurface(bd, 'pial lh');
                bp.loadSurface(bd, 'pial rh');
            end
        end

    end
    
    % Convenience methods
    methods
        function [hems,is_lh,is_rh] = getHems(bd, jack)
            % returns a cell array of hemispheres in the jacksheet (e.g. {'lh'})
            if nargin < 2
                jack = bd.docs.jacksheet;
            end
            hems = sort(unique(jack.whichHemi));
            hems = hems(~cellfun('isempty',hems));
            is_lh = ismember('lh',hems);
            is_rh = ismember('rh',hems);
        end
    end
    
    % Static methods
    methods (Static, Access = public)
        
        function bds = xLoadSubj(subjs, rootEEGdir)
            for i = 1:numel(subjs)
                bds(i) = braindata2(subjs{i}, rootEEGdir);
                bds(i).load_docs();
                bds(i).load_roi();
            end
        end
        function [xDataAvg, minSubjMask, dataMat] = xLeadData2ROI(bds, data, chanNames, varargin)
            % xLeadData2ROI calls leadData2ROI for each subject and averages at each ROIC
            % 
            %
            % INPUT
            %   bds         - array of braindata2 objects
            %   data        - cell array of data per channel
            %   chanNames   - cell array of chan names
            %
            % OPTIONAL KEY-VALUE INPUT
            %   minSubj     - minimum number of subjects to record an ROI. Default 3.
            %
            % OUTPUT
            %   xDataAvg    - average at each ROIC as [lh; rh]
            %   minSubjMask - logical mask of ROIs that meet the minSubj requirment
            %   dataMat     - ROIC x subj matrix of data
            
            ip = inputParser;
            ip.addParameter('rroi',9);
            ip.addParameter('minSubj', 3); 
            ip.parse(varargin{:});
            %rroi    = ip.Results.rroi;
            minSubj = ip.Results.minSubj;
            
            assert(iscell(data) && iscell(chanNames) && numel(data)==numel(chanNames), ...
                'data and chan names must be cells of the same size');
            assert(all(arrayfun(@(x) isa(x, 'braindata2'), bds)) && numel(bds)==numel(data), ...
                'bds must be a braindata2 array of the same size as data');
            
            xDataAvg    = [];
            minSubjMask = []
            
            n = numel(data);
            dataMat = [];
            modargin = [varargin {'force_bilat' 1}];
            for i = 1:n
                [~, valPerROI, sparseMask, rh_begin] = bds(i).leadData2ROI(data{i}, chanNames{i}, modargin{:});
                
                % Initialize
                if i == 1
                    dataMat = nan(length(sparseMask), n);
                end
                
                if rh_begin > 1
                    nroi = length(sparseMask) / 2;
                    dataMat(sparseMask(1:nroi), i)      = valPerROI(1:rh_begin-1);
                    dataMat(sparseMask(nroi+1:end), i)  = valPerROI(rh_begin:end);
                else
                    dataMat(sparseMask, i) = valPerROI;
                end
            end
            
            minSubjMask = (sum(~isnan(dataMat),2) >= minSubj) > 0;

            xDataAvg = nanmean(dataMat, 2);
            
            dbg = 0;
            if dbg
                fprintf('%d lh\n', sum(~isnan(dataMat(1:nroi))));
                fprintf('%d rh\n', sum(~isnan(dataMat(nroi+1:end))));
            end
        end
        function peek(t, n) % Use this for examining your bd.roi.roic_roi table
            % Preview first n (default 10) rows of a table
            
            if nargin < 2, n = 10; end
            if ischar(t), t = evalin('caller',t); end

            n = min(n, height(t));
            disp(t(1:n, :));
            fprintf('\t...%d more rows...\n', height(t) - n);
        end
        function [left, right, lndx, rndx] = chans2hem(chans, jacksheet)
            % chans2hem(chans, jacksheet) takes cell string chans and jacksheet table
            % Input:
            %   chans   - cell string of channel names
            %   jacksheet - the jacktable
            %
            % Output:
            %   left    - left channels 
            %   right   - right channels
            %   lndx    - index into chans
            %   rndx    - index into chans
            
            if any(cellfun(@(x) strfound(x,'-'), chans)) || size(chans,2) == 2
                % for bipolar, it is sufficient to look at the first of the pair
                chans_mp = chans(:,1);
                chans_mp = cellfun(@(x) strsplit(x,'-'), chans_mp, 'uniformOutput',0);
                chans_mp = cellfun(@(x) x{1}, chans_mp, 'uniformOutput',0);
            else
                chans_mp = chans;
            end
            if ismember('PHYS', jacksheet.Properties.VariableNames)
                jacksheet = jacksheet(strcmpi(jacksheet.chanType, 'PHYS'), :);
            end
            jackl = jacksheet(strcmpi(jacksheet.whichHemi, 'lh'), :);
            jackr = jacksheet(strcmpi(jacksheet.whichHemi, 'rh'), :);
            lndx = ismember(chans_mp, jackl.chanName);
            rndx = ismember(chans_mp, jackr.chanName);
            left  = chans(lndx, :);
            right = chans(rndx, :);
        end
        function [subdurals, ndx] = chans2subdural(chans, jacksheet)
            % chans2subdural(chans, jacksheet) takes cell string chans and jacksheet table
            % Input:
            %   chans   - cell string of channel names
            %   jacksheet - the jacktable
            %
            % Output:
            %   subdurals - subdurals channels 
            %   ndx    - index into chans
            
            jack = jacksheet(strcmpi(jacksheet.chanType, 'PHYS'), :);
            [~, depthNdx] = braindata2.chans2depth(chans, jacksheet);
            ndx = ~depthNdx & ismember(chans, jack.chanName);
            subdurals = chans(ndx);
        end
        function [depths, ndx] = chans2depth(chans, jacksheet)
            % chans2depth(chans, jacksheet) takes cell string chans and jacksheet table
            % Input:
            %   chans   - cell string of channel names
            %   jacksheet - the jacktable
            %
            % Output:
            %   depths - depth channels 
            %   ndx    - index into chans
            
            is_2_column_bp = ~isvector(chans) && size(chans,2) == 2;
            assert(isvector(chans) || is_2_column_bp, 'chans must be an n X 1 (or n X 2 for bipolar) cell array');

            
            if is_2_column_bp || any(cellfun(@(x) strfound(x,'-'), chans))
                % for bipolar, it is sufficient to look at the first of the pair
                chans_mp = chans(:,1);
                chans_mp = cellfun(@(x) strsplit(x,'-'), chans_mp, 'uniformOutput',0);
                chans_mp = cellfun(@(x) x{1}, chans_mp, 'uniformOutput',0);
            else
                chans = chans(:);
                chans_mp = chans;
            end

            jack = jacksheet(strcmpi(jacksheet.chanType, 'PHYS'), :);
            jack = jack(strcmpi(jack.hardwareType, 'DEPTH'), :);
            ndx = ismember(chans_mp, jack.chanName);
            depths  = chans(ndx, :);
        end
        function version()
            fprintf('Version: Braindata 2.0 (07/2018)\n');
        end
        
    end
    
    methods (Static) % (Access = private) should technically be private but keeping public for ease of development

        
                
    end
    
    
end
