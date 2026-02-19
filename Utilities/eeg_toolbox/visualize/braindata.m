classdef braindata < handle
% braindata stores and plots a cortical surface and its electrodes.
% The class provides useful operations for viewing and creating figures.
%
% See each of the following for more details:
%
% braindata Properties:
%   meta    - metadata about subject. Meta contains fields: 
%               rootEEGdir
%               hemi
%               subj
%               suma_dir
%
%   surf    - surface struct. Surf contains fields:
%               vertices
%               hemi
%               surfRez
%               faces
%               name
%               surfType
%               subjId
%               dural
%               rmv_medial
%               
%   dural   - dural surface struct
%   curv    - curvature struct
%   parc    - atlas info?
%   leads   - pial xyz coordinates table
%   element_info - element_info table
%   Jacksheet - jacksheetMaster table
%   Monopolar -  monopolar.mat struct
%   FigureData - struct with figure and plot data
%   TagNames - cell array of tag names
%   Hardware - channel names grouped by tag name
%
% braindata Methods:
%   braindata   - constructor
%   plot        - plots a brain
%   plotPoint   - Plots points on the current figure
%   clearPoints - removes all points on current figure
%   loadStandard - Initializes the Colin_n27 standard
%   loadFSAverage - Initializes the fsAverage brain
%   loadSubject - Initializes a subject brain
%   plotHardware - Plots brain and electrodes colored by hardware
%   setOpacity  - Changes surface opacity
%   plotDepths  - Plots only depth hardware
%   plotEcog    - Plots only cortical hardware
%   legend      - Adds legend to figure
%   flip        - Flips brain upside-down
%   plotRegion  - Colors the given vertex set
%   plotData    - Plots 
%   fn_load_surf - loads a surface
%   fn_load_parc - loads the ___
%   fn_whereami_parc - does ____
%   
%
% Author: MST
%
% See Also: braindata2
    methods (Static)
        function version()
            fprintf('Version: Braindata v1.0(01/2018). (See braindata2 for v2.0)\n');
        end
    end
    
    properties (Access = public)
        meta;   % Metabd for the object
        surf;   % A plotable representation of the surface
        dural;
        curv;   % Gives the surface it's shape
        parc;
        leads;  % xyz-coord per lead table (dural)
        element_info;   % element_info table
        Jacksheet;      % jacksheetMaster table
        Monopolar;      % monopolar.mat struct
        FigureData;     % Struct with figure and plot bd
        TagNames;       % cell array of tag names
        Hardware;       % channel names grouped by tag name
        roi;            % struct with fields dists, file, denisity
    end
    
    properties (Access = private)
        subjNum;
        rng;
        
    end
    
    
    
    properties (Constant, Access = private)
        avgSurfDir      = 'data/fsaverage/SUMA'; % relative to braindata
        stdSurfDir      = 'data/colin_N27/suma_MNI';
        stdNumVert      = 198812;                % per hemisphere
        pathSuma        = 'tal/zloc/freesurfer/%s/SUMA';
        pathSumaBackup  = 'tal/intermediates/images_2_fsSumaStd/%s/SUMA';
        pathLeads       = 'tal/leads.csv';
        pathLeadsBackup = 'tal/intermediates/locs_4_compiled/compiled.csv';
        pathMonopolar   = 'tal/monopolar.mat';
        pathRegistered  = 'tal/intermediates/locs_1_registered/coords.mr_pre.ict2bmrAF.csv';
        pathTal         = 'tal';
        pathEdges       = 'tal/intermediates/locs_2_projected/edges.mat';
        pathAnchors     = 'tal/intermediates/locs_0_slicer/coords.mr_pre.anchors.csv';
        pathRoi         = 'tal/intermediates/roi_1_grow';
        ncenter_course  = 600;
        ncenter_fine    = 2400;
    end

    methods
       
        
        function bd = braindata(varargin)
        % braindata accepts the following (optional) key, value pairs:
        %   subj
        %   rootEEGdir
        %   suma_dir
        % These can also be set after initialization using property assignment
        % or with loadSubject
            meta_defaults = {'subj','','parameter';
                'rootEEGdir','','parameter';
                'suma_dir','./','parameter';
                'util2dir',false,'parameter'
                };
            bd.meta = util_extract_inputs(meta_defaults, varargin);
            
            %try bd.element_info = getElementInfo(bd.meta.subj, bd.meta.rootEEGdir); catch; end
            
            try
                bd.rng = rng('shuffle');
            catch e
                rng('default');
                bd.rng = rng('shuffle');
            end

            % add data folder to path
            path_visualize = fileparts(which('braindata'));
            path_data = fullfile(path_visualize, 'data');
            pathCell = regexp(path, pathsep, 'split');
            if ~any(strcmpi(path_data, pathCell))
                addpath(genpath(path_data));
            end
        end
        
        function loadSubject(bd, subj, rootEEGdir, varargin)
        % LOADSUBJECT loads patient specific bd, including the surface
        %   Optional input
        %     subj
        %     rootEEGdir
        %   Optional key-value pairs
        %     hemi, (rh/lh)
        %     surfType, (pial,inflated,inf_200,smoothwm,sphere,sphere.reg,white)
        %     rmv_medial, (true,[false])
        %
        %
        % This should be called whenever we're plotting a subject brain
        % Important: if your subject has electrodes on both hemispheres,
        % you should only load in one. Use the 'hemi','lh'/'rh' param pair to
        % do this. In the future, this class will work without you needing to do
        % this.
        %
            bd.meta.subj = subj;
            bd.meta.rootEEGdir = rootEEGdir;
            
            %fprintf('Subject %s\n', subj);
            [~,bd.meta.subjNum] = util_split_stringnum(subj);
            if isnan(bd.meta.subjNum), bd.meta.subjNum = randi(100,1); end
            bd.rng = rng(bd.meta.subjNum); %#ok<*CPROPLC>
            
            % define default inputs
            varargin = util_ensure_varargin(varargin);
            default_cell = {'hemi',[],'parameter';
                            'surfType','pial','parameter';
                            'rmv_medial', 0, 'parameter'};
            inputs = util_extract_inputs(default_cell, varargin);
            hemi = inputs.hemi;
            
            % set suma path
            fname = fullfile(rootEEGdir, subj, bd.pathSuma);
            if ~isempty(strfind(fname, '%s')), fname = sprintf(fname, subj); end
            if ~exist(fname, 'dir')
                temp = sprintf(bd.pathSumaBackup, subj);
                fname = fullfile(rootEEGdir, subj, temp);
            end
            if exist(fname, 'dir')
                bd.meta.suma_dir = fname;
                %fprintf('Suma dir set\n');
            end
        
            % load element info
            try 
                bd.element_info = getElementInfo(bd.meta.subj, bd.meta.rootEEGdir);
                bd.TagNames = bd.element_info.tagName;
                bd.TagNames = bd.TagNames(ismember(bd.TagNames, util_split_stringnum(getLeads(bd.meta.subj, bd.meta.rootEEGdir))));
                %fprintf('Element info loaded\n');
            catch e
                fprintf('Error loading element info: %s\n', e.message);
            end
            
            % load jacksheet
            try 
                bd.Jacksheet = getJackTable(bd.meta.subj, bd.meta.rootEEGdir); 
                %fprintf('jacksheet loaded\n');
            catch e
                fprintf('Error loading jacksheet %s\n', e.message);
            end
            
            % load leads xyz
            fname = fullfile(rootEEGdir, subj, bd.pathLeads);
            if ~exist(fname, 'file')
                fname = fullfile(rootEEGdir, subj, bd.pathLeadsBackup);
            end
            if exist(fname, 'file')
                bd.leads = readtableSafe(fname);
                leadsGot = getLeads(bd.meta.subj, bd.meta.rootEEGdir);
                if ~isempty(leadsGot)
                    bd.leads = bd.leads(ismember(bd.leads.chanName, leadsGot), :);
                end
                leads = bd.leads;
                if ~isempty(bd.TagNames) && ~isempty(leads) && ~isempty(bd.element_info)
                    bd.Hardware = groupBySubstring(leads.chanName, bd.TagNames);
                    bd.FigureData.colors = rand(length(bd.Hardware), 3);
                end
                %fprintf('Leads loaded (%s)\n', fname);
            else
                fname = fullfile(rootEEGdir, subj, bd.pathRegistered);
                if exist(fname, 'file')
                    leads = readtableSafe(fname);
                    if ~isempty(bd.TagNames) && ~isempty(leads) && ~isempty(bd.element_info)
                        bd.Hardware = groupBySubstring(leads.chanName, bd.TagNames);
                        bd.FigureData.colors = rand(length(bd.Hardware), 3);
                    end
                    %fprintf('Leads loaded (%s)\n', fname);
                end
            end
            
            % load monopolar
            fname = fullfile(rootEEGdir, subj, bd.pathMonopolar);
            if exist(fname, 'file')
                mbd = load(fname);
                bd.Monopolar = mbd.monopolar;
                %fprintf('monopolar loaded\n');
            end
            
            % load surface
            bd.loadSurf(hemi, inputs.surfType, inputs.rmv_medial);
            
        end
        
        function loadSurf(bd, hemi, surfType, rmv_medial, suma_dir)
            if ~exist('rmv_medial','var') || isempty(rmv_medial), rmv_medial = 0; end
            if exist('suma_dir', 'var')
                
            end
            
            if isempty(hemi)
                % If subj has electrodes on only 1 hemisphere, load that one
                if ~isempty(bd.element_info)
                    hemis = unique(bd.element_info.whichHemi);
                    hemis = hemis(~cellfun(@isempty,hemis));
                else
                    hemis = {'rh','lh'};
                end

                bd.meta.hemi = hemis;
                
                if length(hemis) == 1
                    bd.fn_load_surf('hemi',char(hemis), 'surfType',surfType, 'rmv_medial',rmv_medial);
                else
                    % otherwise load both
                    bd.fn_load_surf('surfType',surfType,'rmv_medial',rmv_medial);
                end
            else
                % use specified a hemisphere
                bd.fn_load_surf('hemi', hemi, 'surfType',surfType,'rmv_medial',rmv_medial);
                bd.meta.hemi = hemi;
            end
        end
        
        function parcellate(bd)
        % PARCELLATE colors the brain based on Destreaux atlas labels
            if isempty(bd.parc)
                d = bd.meta.suma_dir;
                hemi = bd.meta.hemi;
                args = {'dir',d};
                if ~isempty(hemi) && (ischar(hemi) || (iscellstr(hemi) && isscalar(hemi)))
                    args = [args {'hemi',char(hemi)}];
                end
                bd.fn_load_parc(args{:});
            end
            
            try surf = findall(bd.FigureData.figure,'type','patch'); surf=surf(1);
            catch e, surf=[]; end
            
            lut = bd.parc.a2009s.lut;
            labels = lut.label_name;
            rgb = cat(2, lut.r, lut.g, lut.b);
            verts = lut.vert_idx; % cells of vertex arrays
            bd.surf.rgb = nan(size(bd.surf.vertices,1), 3);
            
            
            for i = 1 : length(verts)
                vi = verts{i};
                rgbi = double(repmat(rgb(i,:), length(vi), 1)) / 255;
                bd.surf.rgb(vi,:) = rgbi;
                if ~isempty(surf)
                    if size(surf.FaceVertexCData,3) ~= 3
                        surf.FaceVertexCData = bd.surf.rgb;
                    end
                    surf.FaceVertexCData(vi,:) = rgbi;
                end
                
                displayName = atlasFSTranslate('destrieux','surf',labels{i});
                if isempty(displayName)
                    x=patch(nan,nan,rgbi(1,:));
                else
                    x=patch(nan,nan,rgbi(1,:),'displayName',displayName);
                end
            end
            
            
            % params for vutil_prep_layers
            % plot() defines only ulay_dset
            
%             'num_pts',   [], 'parameter';
%         'olay_rgb',  [], 'parameter';
%         'ulay_rgb',  [], 'parameter';
%         'olay_dset', [], 'parameter';
%         'ulay_dset', [], 'parameter';
%         'olay_cmap', struct(), 'parameter';
%         'ulay_cmap', struct(), 'parameter';
%         'olay_thresh',[NaN NaN],'parameter';
%         'blend', 0, 'parameter'};
            
        end
 
        function fn_load_surf(bd, varargin)
        % FN_LOAD_SURF needs to be called to load in the surface
        %
        % Designed to work with SUMA/ files
        % Expected filenames are like: lh.pial.gii or std.141.rh.pial.gii
        %  which comes approx from [surfRez].[hem].[surfType].gii
        %
        % The following (optional) key, value pairs are accepted*:
        %   subj                      - will use class property if not passed
        %   hemi {[separate], combined, lh, rh} - load a single hemisphere.
        %                                 Combined loads both, separate returns a
        %                                 struct with an lh and rh braindata properties
        %   surfRez {[std.141], native} - how many vertices to expect. This also
        %                                 determins what filename to look for
        %   surfType {[pial], inflated} - determines filename
        %   dir                         - where to look for the surface files (usually the SUMA dir)
        %   rmv_medial                  - ? (Mike doesn't know what this does yet)
        %   dural {true, [false]}       - whether to load the dural surface
        %                                 (ie pial-outer-smoothed) in addition to the pial
        %
        % *{options} [default]
            bd.surf = []; % reset bd.surf
            bd.curv = [];
            
            % define default inputs
            varargin = util_ensure_varargin(varargin);
            default_cell = {...
                'subj',bd.meta.subj,'parameter';
                'hemi','separate','parameter';
                'surfRez','std.141','parameter';
                'surfType','pial','parameter';
                'util2dir',bd.meta.util2dir,'parameter';
                'dir',bd.meta.suma_dir,'parameter';
                'rmv_medial',false,'parameter';
                'dural',false,'parameter'
                };
            inputs = util_extract_inputs(default_cell, varargin);
            %hemi = inputs.hemi;
            %bd.surfRez = inputs.surfRez;
            fprintf('Grabbing surface as %s hemisphere(s)...\n', inputs.hemi);
            
            
            % load the actual surface
            %fprintf('\tLoading surface...\n');
            inputs.rootEEGdir = bd.meta.rootEEGdir;
            bd.surf = surf_load_main(inputs);

            % remove imaginary medial surface, if so desired
            if inputs.rmv_medial
                fprintf('\tRemoving imaginary medial surface...\n');
                bd.fn_load_parc(inputs)
                sub_surf = surf_extract_parcel(bd.surf, ...
                    bd.parc.a2009s.lut.label_name, ...
                    bd.parc.a2009s.lut.vert_idx, ...
                    {'wm_lh_Unknown', 'wm_rh_Unknown', 'Unknown'}, ...
                    'inverse',1);
                bd.surf.faces = sub_surf.faces;
                bd.surf.num_verts_shown = length(unique(sub_surf.faces(:)));
                
                % Edit: we actually want to keep vertices in place so we don't have to
                % change anything else
                % remove vertices not used in faces
%                 sub_v_ndx = sort(unique(sub_surf.faces(:))); % list 1:n, but with some missings (n-n' missing)
%                 bd.surf.vertices = bd.surf.vertices(sub_v_ndx, :);         % length 1:n'
%                 
%                 % Now we need to replace face indices with new correspondence
%                 bd.surf.faces = indexOf(sub_surf.faces, sub_v_ndx); % map vals 1:n to 1:n'
%                 
                
            end
            
            % load dural 
            if inputs.dural
                fprintf('\tLoading the dural surface...\n');
                dural_inputs = inputs;
                dural_inputs.surfRez = 'native';
                dural_inputs.surfType = 'pial-outer-smoothed';
                bd.dural = surf_load_main(dural_inputs);
            end
            
            % load curvature
            try
                %fprintf('\tLoading curvature...\n');
                inputs.dset = 'sulc';
                bd.curv = surf_load_dset(inputs);
            catch
                fprintf('\t\tCurvature not found!\n')
            end
            
            %fprintf('Surface loaded!\n');
        end
        
        function loadFSAverage(bd, varargin)
        % LOADFSAVERAGE loads the fsaverage surface into braindata
        %   This function will call fn_load_surf with the appropriate 'dir'
        %   and 'subj' arguments set.
        %
        %   *Any additional arguments you give it will be passed to fn_load_surf
        %   you are responsible for making sure those params make sense
        %
        %   See also: FN_LOAD_SURF
            
            %
            %bd.meta.subj = 'fsaverage';
            bdDir = fileparts(which('braindata'));
            d = fullfile(bdDir, bd.avgSurfDir);
            newArgs = {'dir', d, 'subj', 'fsaverage'};
            allArgs = [newArgs, varargin];
            bd.fn_load_surf(allArgs{:})
            if isfield(bd.surf, 'hemi')
                bd.meta.hemi = bd.surf.hemi;
            end
        end
        
        function loadStandard(bd, varargin)
        % LOADSTANDARD loads the Colin n27 surface into braindata
        %   This function will call fn_load_surf with the appropriate 'dir'
        %   and 'subj' arguments set.
        %   Any additional arguments you give it will be passed to fn_load_surf
        %   YOU are responsible for making sure those params make sense
        %   See also: FN_LOAD_SURF
            
            %
            bd.meta.subj = 'colin_n27';
            bdDir = fileparts(which('braindata'));
            d = fullfile(bdDir, bd.stdSurfDir);
            newArgs = {'dir', d, 'subj', 'standard_colin_n27'};
            allArgs = [newArgs, varargin];
            bd.fn_load_surf(allArgs{:})
        end
        
        function fn_load_parc(bd, varargin)
        % FN_LOAD_PARC does ____
        % must call with surface loaded as combined
        %
        % Accepts the following (optional) key-value pairs:
        %   subj - will use class property if not passed
        %   hemi - uses surface's property if not overriden
        %   surfRez - uses surface's property if not overriden
        %   parcType [a2009s]
        %   dir
            try
                default_cell = {
                    'subj',bd.meta.subj,'parameter';
                    'hemi', bd.surf.hemi, 'parameter';
                    'surfRez',bd.surf.surfRez,'parameter';
                    'parcType','a2009s','parameter';
                    'util2dir',bd.meta.util2dir,'parameter';
                    'dir',bd.meta.suma_dir,'parameter'
                    };
            catch e
                warning('You may have tried to call load_parc without using a combined surface...');
                fprintf('Caught this error: %s\n', e.message);
                fprintf('I will relad the surface with combined and try it again for you.\n');
                keyboard;
                fprintf('Calling fn_load_surface("hemi","combined..."\n');
                bd.fn_load_surf('hemi','combined');
                bd.fn_load_parc(varargin{:});
                return;
            end
            inputs = util_extract_inputs(default_cell, varargin);
            
            % now load the parc file
            bd.parc.(inputs.parcType) = surf_load_parc(inputs);
        end
        
        function parcs = fn_whereami_parc(bd, nodes, varargin)
        % Written by Cocjin
            iscell_nodes = iscell(nodes);
            if ~iscell_nodes
                nodes = {nodes};
            end
            
            % define default inputs
            default_cell = {
                'parcType','a2009s','parameter';
                'labelField','label_name','parameter'};
            inputs = util_extract_inputs(default_cell, varargin);
            
            % parcs per thing.
            parcs = surf_nodes2parcels(bd.parc.(inputs.parcType).lut.(inputs.labelField),...
                bd.parc.(inputs.parcType).lut.vert_idx,...
                nodes);
            
            if ~iscell_nodes
                parcs = parcs{1};
            end
        end
        
        function overlayDura(bd)
            if isempty(bd.dural) && ~isempty(bd.surf)
                bd.fn_load_surf('hemi',bd.surf.hemi,'dural',1);
            end
            h = bd.plot('dural',1,'handle',bd.FigureData.figure);
        end
        
        function h = plot(bd, varargin)
        % plots a brain bd object using viz_surf_main
        %
        % Optional (key-value) Inputs:
        %   name - Provide a name for the figure window. default: meta.subj
        %   snapCurv - true/[false]. Snaps negative/positive curvatures to
        %               -1/1. Hint: pass camlight=true for better shading.
        %   view - One of these strings: left-sagittal, right-sagittal,
        %          dorsal, ventral. Default is left-sagittal.
        %   hemi - Either lh or rh
        %   handle - pass in a figure handle to plot onto. (otherwise new figure)
        %   newAxis - True/[False]. True to create plot with new axis
        %   camlight - True/[False]. True to add camlight
        %   dural - True/[False]. Whether to load the dural
        %   surfType - If you want to plot a surface that's not currently loaded. See surfTypes.
        %   skipCurv - True/[False] if you want to skip curvature

        % Outputs
        %   h - figure handle
        
            % argument handling
            ip = inputParser;
            ip.addParameter('name', bd.meta.subj, @ischar);
            ip.addParameter('snapCurv', false);
            ip.addParameter('view', [], @ischar);
            ip.addParameter('hemi', '', @ischar);
            ip.addParameter('handle', []);
            ip.addParameter('newAxis', false, @islogical);
            ip.addParameter('camlight',false);
            ip.addParameter('dural',false);
            ip.addParameter('surfType',[]);
            ip.addParameter('skipCurv',0);
            ip.KeepUnmatched = true;
            ip.parse(varargin{:});
            name = ip.Results.name;
            snapCurv = ip.Results.snapCurv;
            camView = ip.Results.view;
            h = ip.Results.handle;
            hemi = ip.Results.hemi;
            newAxis = ip.Results.newAxis;
            useLight = ip.Results.camlight;
            dural = ip.Results.dural;
            surfType = ip.Results.surfType;
            skipCurv = ip.Results.skipCurv;
            
            crv = [];
            
            % preconditions
            assert(~isempty(bd.surf), 'You need to load the surface first');
            assert(~isempty(bd.curv) || skipCurv, 'No curvature - Try loading the surface first');
            
           % Possibly load specified surface
            if ~isempty(surfType) && ~strcmpi(surfType, bd.surf.surfType)
                bd.loadSurf(hemi, surfType);
            end
            
            % setup
            surf = bd.surf;
            if isempty(hemi) && isempty(camView) && isfield(bd.meta, 'hemi')
                if strcmpi(bd.meta.hemi, 'lh')
                    camView = 'left-sagittal';
                else
                    camView = 'right-sagittal';
                end
            elseif isempty(camView)
                camView = 'right-sagittal';
            end
            
            if isempty(h)
                h = figure();
                newFigure = true;
            else
                newFigure = false;
            end
            fieldsSurf = fieldnames(bd.surf);
            if ismember('hemi', fieldsSurf)
                if ~isempty(hemi) && ~strcmpi(bd.surf.hemi, hemi)
                    error('Surf is loaded without hemisphere you specified');
                end
                hemi = bd.surf.hemi;  %#ok<*NASGU>
                if dural
                    surf = bd.dural;
                elseif ~skipCurv
                    surf = bd.surf;
                    crv = bd.curv.data;
                end
                
                isDuralCurv = ~skipCurv && bd.curv.dural;
                
            elseif ismember('lh', fieldsSurf) && ismember('rh', fieldsSurf)
                if ~isempty(hemi)
                    if dural
                        surf = bd.dural.(hemi);
                    elseif ~skipCurv
                        surf = bd.surf.(hemi);
                        crv = bd.curv.(hemi).data;
                    end
                    
                    isDuralCurv = bd.curv.(hemi).dural;
                else
                    bd.plot(varargin{:}, 'handle', h, 'hemi', 'lh');
                    bd.plot(varargin{:}, 'handle', h, 'hemi', 'rh', 'newAxis',false);
                    return;
                end
            end
            
            if dural && ~isDuralCurv
                warning('The curvature loaded is not dural; you need to load the dural first!');
            end
            
            %crv = bd.curv.data;
            if snapCurv && ~skipCurv
                crv(crv < 0) = -1;
                crv(crv > 0) = 1;
                % take this following line out!
                %crv = []; %ones(length(crv),1);
            end
            
            
            
            % figure creation
            if ~newFigure && newAxis
                % add new axis
                axes;
            end
            h.Name = name;
            if ismember('rgb',fieldnames(bd.surf)) && ~isempty(bd.surf.rgb)
                olay_rgb = bd.surf.rgb;
                %surfPatch = viz_surf_main(surf, 'ulay_dset', crv, 'olay_rgb',olay_rgb);
                surfPatch = viz_surf_main(surf, 'ulay_dset', crv);
            else
                surfPatch = viz_surf_main(surf, 'ulay_dset', crv);
            end
            viz_set_shine(surfPatch);
            axes(surfPatch.Parent);
            axis equal;
            axis tight;
            shading interp;
            lighting gouraud;
            material dull;
            axis off;
            hold on;
            
            if dural
                surfPatch.FaceAlpha = 1;
                surfPatch.FaceColor = [0 0 1];
            end
            
            %zoom(3);
            switch lower(camView)
                case 'ventral'
                    view(180,-90);
                case 'left-sagittal'
                    view(-90,0);
                case 'right-sagittal'
                    view(90,0);
                case 'dorsal'
                    view(0,90);
            end
            
            bd.FigureData.figure = h;
            if useLight
                camlight infinite
            end
        end
        
        function surfTypes(bd)
        % Display all surf types that can be loaded to this braindata
            d1 = fullfile(bd.meta.rootEEGdir, bd.meta.subj, sprintf(bd.pathSuma,bd.meta.subj));
            d2 = fullfile(bd.meta.rootEEGdir, bd.meta.subj, sprintf(bd.pathSumaBackup,bd.meta.subj));
            f1 = lsCell(sprintf('%s/std.141.*.gii', d1));
            f2 = lsCell(sprintf('%s/std.141.*.gii', d2));
            files = union(f1, f2);
            [~,files,types] = cellfun(@fileparts, files, 'uniformoutput',0);
            [~,~,types] = cellfun(@fileparts, files, 'uniformoutput',0);
            types = cellfun(@(x) x(2:end), types, 'uniformoutput',0);
            types = unique(types);
            disp(types);
            
        end
        
        function plot2(bd)
            bd.plot; bd.setOpacity(0.6); bd.clearPoints;
            % color by tag
            tags = bd.TagNames;
            tags = tags(~cellfun('isempty',tags));
            chans = bd.leads.chanName;
            xyz = bd.leads(:,{'x','y','z'});
            [~,ndxs] = groupBySubstring(chans, tags);
            colors = hsv(numel(tags));
            for i = 1:numel(ndxs)
                chan = chans(ndxs{i});
                [~,num] = util_split_stringnum(chan);
                bd.plotPoint(xyz(ndxs{i},:), 'color', colors(i,:), 'legend', tags{i}, 'label',arrayfun(@(x) {x}, num));
            end
            bd.legend();
            figfmt;
        end
        
        function [hgroup, hsphere] = plotPoint(bd, points, varargin)
        % Plots points on the current figure
        %
        % If plotting a subject's point on their own surface, you can pass XYZ
        % coordinates. If plotting onto an avergae or other standard surface, pass
        % indices.
        %
        % Inputs:
        %   points - Either a list of xyz coords, a list of mesh indices, or a list
        %           of channelNames, or a table with x,y,z columns
        % Optional (key-value) Inputs:
        %   color - default [0 1 0] green for each point
        %   radius - default 1.5 (mm)
        %   legend - string label that applies to all the points passed.
        %           Sets the display name of the group 
        %   ctNum - If you pass in 1 or 2, when looking up xyz-coords for channel
        %           names, will look in CT_[ctNum]/leads.csv instead of
        %           tal/leads.csv
        %   useNode - default false. If true, plot at vertex instead of xyz
        %   label - Plot the text labels given in a cell array
        %   alpha - alpha level
        %   data  - mapped into colormap
        %   hemi - If surf is a struct with lh and rh fields, pass 'lh'/'rh' to determine plot
        %
        % Outputs:
        %   hgroup - Handle to group of the point sphere(s)
        %   hsphere - Handle(s) to point sphere(s)
        %
        %   
        %
        % Note that if you pass mesh indices, points will be on the pia.
        % If you pass xyz's that coorespond to dural points, the points may be just
        % above the brain's pial surface
        %
        % If you pass channel names, this function will look in tal/leads.csv for
        % coordinates
        %
        
            if istable(points) 
                points = cat(2,points.x, points.y, points.z);
            end
            
            n = length(points);
            if n==0, return; end
            subj = bd.meta.subj;
            rootEEGdir = bd.meta.rootEEGdir;
            
            
            % argument handling
            ip = inputParser;
            ip.addParameter('color', repmat([0 1 0], n, 1));
            ip.addParameter('radius', 1.5, @isnumeric);
            ip.addParameter('legend', '', @ischar);
            ip.addParameter('ctNum', [], @isnumeric);
            ip.addParameter('useNode', false);
            ip.addParameter('label', '');
            ip.addParameter('alpha',1);
            ip.addParameter('data',[]);
            ip.addParameter('hemi',[]);
            ip.KeepUnmatched = true;
            ip.parse(varargin{:});
            color = ip.Results.color;
            radius = ip.Results.radius;
            legendLabel = char(ip.Results.legend);
            ctNum = ip.Results.ctNum;
            useNode = ip.Results.useNode;
            label = ip.Results.label;
            alpha = ip.Results.alpha;
            data = ip.Results.data;
            hemi = ip.Results.hemi;
            
            if ischar(points), points = cellstr(points); end
            arePointsNames = iscellstr(points);
            
            if isrow(points) && numel(points) ~= 3
                points = points';
            end
            
            % preconditions
            assert(~isempty(bd.surf), 'You need to load the surface first');
            if ~arePointsNames
                assert(size(points,2) == 1 || size(points,2) == 3, 'plotPoint takes points with size n-by-1 or n-by-3'); 
            end
            
            
            
            
            % convert indices to xyz, and tagNames to xyz
            if arePointsNames
                    points = cellstr(points);
                    
                    % Use nodes passed or should be done b/c surf isnt pial
                    if isempty(strfind(bd.getASurf().name, 'pial')) || useNode
                        fprintf('Using a std vertex for point\n');
                        nodes = bd.getChanNodes(points,1);
                        bd.plotPoint(nodes, varargin{:});
                        return;
                        
                    end
                    
                    % load leads table
                    if isempty(ctNum)
                        leadsTable = bd.leads;
                    else
                        filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, sprintf('tal/CT_%d/leads.csv',ctNum));
                        if ~exist(filename, 'file')
                            fprintf('No leads file found: %s\n', filename);
                            return
                        end
                        leadsTable = readtableSafe(filename);
                    end
                    
                    % leave this code in case we switch to preferring monopolar
%                     if ~isempty(bd.Monopolar)
%                         rowMask = ismember({bd.Monopolar.chanName}, points);
%                         if isempty(find(rowMask, 1)) 
%                             fprintf('No channels found\n');
%                             return;
%                         end
%                         temp = {bd.Monopolar(rowMask).native_xyz};
%                         temp = cellfun(@(t) [t.x t.y t.z], temp, 'uniformOutput',false);
%                         xyz = cell2mat(temp');
                    if ~isempty(leadsTable)
                        rowMask = ismember(leadsTable.chanName, points);
                        if isempty(find(rowMask, 1)) 
                            fprintf('No channels found (hint: did you forget an "s" for shift?\n');
                            return;
                        end
                        rows = leadsTable(rowMask,:);
                        if ~isempty(label)
                            label = label(sortBySubstring(label, leadsTable.chanName));
                        end
                        xyz = [rows.x, rows.y, rows.z];
                        
                    elseif ~isempty(subj) && ~isempty(rootEEGdir)
                        xyz = [];
                        for i = 1 : length(points)
                            loc = getLocation(subj, rootEEGdir, points{i});
                            xyz = [xyz; loc]; %#ok<AGROW>
                        end
                    else
                        fprintf(['Looks like you passed me channel names but I'...
                            'dont know where to look for them - try loading a subject\n']);
                    end
            else
                if size(points, 2) == 1
                    % handle vertices
                    if ~ismember('vertices', fieldnames(bd.surf))
                        % must have 'lh' and 'rh' fields
                        if isempty(hemi)
                            % this is handling for if you are using a convention of passing twice the index range
                            % ie the 198,813'th point is the 1st point on the right hemisphere
                            lh_points = points(points <= bd.stdNumVert);
                            rh_points = points(points > bd.stdNumVert);
                            xyz = [bd.surf.lh.vertices(lh_points, :); bd.surf.rh.vertices(rh_points)];
                        else
                            xyz = bd.surf.(hemi).vertices(points,:);
                        end
                    else
                        xyz = bd.surf.vertices(points, :);
                    end
                    
                else
                    % size is 3 => must be coords
                    xyz = points;
                end
            end
            
            
            
            % plot
            [hgroup, hsphere] = viz_view_sphere(xyz, 'color', color, 'fill', color, 'radius', radius, 'figure', bd.getFigure, 'opacity',alpha);
            if ~isempty(data) && numel(data) == 1
                for k = 1:length(hsphere)
                    hsphere.CData = repmat(data, size(hsphere.XData));
                    hsphere.FaceColor = 'flat';
                end
            end
            
            if ischar(label) && ~iscellstr(label), label = {label}; end
            for i = 1 : length(label)
                cam = get(gca, 'CameraPosition');
                v = normalise(cam - xyz(i,:)); % vector from point to cam
                label_pos = xyz(i,:) + 2*v; % this pops label out a bit so it's visible
                t = text('Position',label_pos, 'String',label{i}, 'FontSize',13, 'HorizontalAlignment','center');
                t.FontWeight = 'bold'; 
            end
            
            if ~isempty(legendLabel)
                hgroup.DisplayName = char(legendLabel);
            end
        end
        
        function f = getFigure(bd)
            if ~isempty(bd.FigureData) && isfield(bd.FigureData,'figure')
                f = bd.FigureData.figure;
            else
                f = gcf;
            end
        end
        
        function index = xyzToVertex(bd, xyz, sphereModel)
        % XYZTOVERTEX findes the nearest standard 141 index to the given coordinate
        %
        % This should ONLY be used to plot Utahs. Utahs are modeled as single
        % points. Other electrodes are modeled as a sphere, so there are
        % usually multiple mesh nodes that correspond to them.
        % 
        % Inputs:
        %   xyz - xyz coordinates (may pass in multiple)
        %
            if nargin < 3, sphereModel = 0; end
            
            if sphereModel
                dural_xyz = [];
                pail_xyz = [];
                r = [];
                snap_to_dura_first = 0;
                indices = lead_to_pial([]);
                index = indices{1};
            else
                assert(~isempty(bd.surf), 'You need to load the surface first');
                n = size(xyz, 1);
                index = zeros(n, 1);
                for i = 1 : n
                    [ndx, dist] = knnsearch(bd.surf.vertices, xyz(i, :), 'k',1);
                    fprintf('Point %d: vertex found %d mm away\n', i, dist(1));
                    index(i) = ndx(1);
                end
            end
            
        end
        
        function nodes = chansToVertex(bd, chanNames, single)
        % chanToVertex findes the nearest standard 141 index to the given channel
        %
        % 
        % Inputs:
        %   xyz - xyz coordinates (may pass in multiple)
        %   single - optional, default=0, if 1, give only 1 node, if 0 give all
        
            if nargin < 3, single = 0; end
            
            nodes = getChanNodes(bd, chanNames, single);
            
            if ~iscellstr(chanNames)
                nodes = nodes{1};
            end

            
        end
        
        function index = reflectPoint(bd, ndxs)
        % REFLECTPOINT projects all points from either hemisphere to one hemisphere
        %
        % This function is useful ONLY for plotting on the fsaverage brain, since
        % the average brain is symmetric. For instance, if we have a set of points
        % from both hemispheres, we may want to show their standardized location on
        % a single figure with a single hemisphere. This will take all indices
        % which are on the other hemisphere and place them in their corresponding
        % location on the contralateral hemisphere using the simple operation of
        % reflecting the X-coordinate.
        %
        % Inputs:
        %   ndxs - Indices on the standard-141 brain
        %
        % Returns index on opposite hemisphere
            
            % preconditions
            assert(~isempty(bd.surf), 'You need to load the surface first');
            assert(0 == mod(length(bd.surf.vertices), bd.stdNumVert), 'Looks like you dont have a std surface');
            
            std_xyz_in = bd.surf.vertices(ndxs, :);   % get xyz
            xyz = [-std_xyz_in(1), std_xyz_in(2:3)];    % flip x
            index = bd.xyzToVertex(xyz);              % get nearest vertex
            
        end
        
        function compare_implant(bd)
            
            d = fileparts(fullfile(bd.meta.rootEEGdir, bd.meta.subj, bd.pathLeadsBackup));
            files = lsCell(d);
            t_leads = [];
            for i = 1:length(files)
                t_leads = [t_leads {readtable(fullfile(d, files{i}))}];
            end
            hemis = unique(bd.Jacksheet.hemi);
            for ih = 1:length(hemis)
                hemi = hemis{ih};
                bd.loadSubject(bd.meta.subj, bd.meta.rootEEGdir, 'hemi', hemi);
                bd.plot();
                for i = 1 : length(files)
                    t = t_leads{i};
                    c = [0 0 0];
                    c(i) = 1;
                    fprintf('\nCT %s color:', files{i});
                    disp(c);
                    for j = 1:height(t)
                        ji = find(strcmpi(bd.Jacksheet.chanName, t.chanName{j}),1);
                        if ji && strcmpi(bd.Jacksheet.hemi{ji}, hemi)
                            bd.plotPoint(t(j,:), 'color', c, 'label',[t.chanName{j}, '_', num2str(i)]);
                        end
                    end
                end
            end
        end
        
        function plotHardware(bd, hw, varargin)
            % plotHardware plots a hardware set on the brain
            % plotHardware(bd, hw, varargin)
            %
            % IMPORTANT: This function is designed to be given a subset of the cells in
            % the braindata property braindata.Hardware
            %
            % Example: bd.plotHardware({bd.Hardware{1}, bd.Hardware{3}});
            % (note bd.Hardware is ordered by bd.TagNames)
            %
            % Inputs:
            %   hw - cell array of chanNames on the hardware piece
            %
            % Optional Input:
            %   Takes all same inputs as plot
            %   registered - use CT-registered coords
            %   ctnum      - numeric, use CT-1 or CT-2
            %   color      - if empty use default
            %
            % Reads from tal/leads.csv
            %
            % See also: plotPoint, plotEcog

            ip = inputParser;
            ip.KeepUnmatched = 1;
            ip.addParameter('registered',0);
            ip.addParameter('ctnum',1);
            ip.addParameter('color',[]);
            ip.parse(varargin{:});
            registered = ip.Results.registered;
            ctnum = ip.Results.ctnum;
            color = ip.Results.color;
            
            % not using pial surf
            if isempty(strfind(bd.getASurf().name, 'pial'))
                fprintf('pial surf not loaded, using std nodes instead of xyz\n');
                useNodes = true;
            else
                useNodes = false;
            end
            
            % if got tagnames, map tagnames to cells of chanNames
            if ischar(hw) || iscellstr(hw)
                if ~iscell(hw), hw = {hw}; end
                hw_tags = hw;
                hw = cell(size(hw_tags));
                lookup_hw = @(tag) bd.Hardware(find(strcmpi(bd.TagNames, tag), 1));
                hw = cellfun(lookup_hw, hw_tags);
            end
            
            % if only 1 hardware piece, ensure it's in a cell
            if ~isempty(hw)
                elmt = hw{1};
                if ~iscell(elmt) && ischar(elmt)
                    fprintf('Expected cell array of cells of chars but got cell array of chars. ');
                    fprintf('Putting your hardware into a cell for you.\n');
                    hw = {hw};
                end
            end
            
            d = fullfile(bd.meta.rootEEGdir, bd.meta.subj, bd.pathTal);
            if registered
                d = fullfile(fileparts(d), bd.pathRegistered);
                if ctnum > 1
                    d = strrep(d, 'mr_pre', ['mr_pre_' num2str(ctnum)]);
                end
            elseif ctnum > 1
                d = fullfile(d, ['CT_' num2str(ctnum)]);
            end
            
            
            
            csv = lsCell(fullfile(d, '*.csv'));
            if length(csv) > 1
                disp(csv);
                choice = str2double(input('Which one? (enter num) ','s'));
                csv = csv(choice);
            end
            
            t = readtableSafe(fullfile(d, char(csv)));

            % get hardware tagNames
            if isempty(hw)
                warning('No hardware on this hemisphere'); 
                bd.plotPoint(t, varargin{:});
                return;
            end
            hw = hw(~cellfun('isempty',hw));
            temp = cellfun(@util_split_stringnum, hw, 'uniformOutput',false); % c array of c array
            tagNames = cellfun(@(x) x{1}, temp, 'uniformOutput',false);
            mask = ismember(bd.TagNames, tagNames);
            colors = bd.FigureData.colors(mask,:);
            assert(length(colors) >= length(hw));
            
            % new plot
            %if ~any(cellfun(@(x) strcmpi(x, 'handle'), varargin))
                bd.plot(varargin{:});
            %end

            chanNames = t.chanName;
            xyz = [t.x, t.y, t.z];
            
            for i = 1 : length(hw)
                hardwareNames = hw{i};
                mask = ismember(chanNames, hardwareNames);
                
                hardware_xyz = xyz(mask,:);
                if isempty(hardware_xyz), continue; end
                lgd = unique(util_split_stringnum(hardwareNames));
                if isempty(color)
                    if isempty(colors)
                        col = repmat(rand(1, 3), size(hardware_xyz,1), 1);
                    else
                        col = colors(i,:);
                    end
                end
                if useNodes
                    bd.plotPoint(chanNames(mask), 'color', col, 'legend', char(lgd), varargin{:});
                else
                    bd.plotPoint(hardware_xyz, 'color', col, 'legend', char(lgd), varargin{:});
                end
            end
            
        end
    
        function plotEcog(bd, varargin)
            % Plots all strips and grids (location type=cortical)
            %
            % Optional Input:
            %   same as plot()
            %
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.addParameter('hemi','',@ischar);

            ip.parse(varargin{:});
            hemi = ip.Results.hemi;
            
            if isempty(hemi) && isfield(bd.meta, 'hemi') && (ischar(bd.meta.hemi) || length(bd.meta.hemi) == 1)
                hemi = bd.meta.hemi;
            end
            
            
%             hemiMask = ones(length(bd.Hardware), 1);
%             if ~isempty(hemi)
%                 assert(ismember(hemi, {'lh','rh'}));
%                 hemiMask = strcmpi(bd.element_info.whichHemi, hemi);
%             end
%             
%             
%             if isempty(bd.element_info)
%                 locMask = ones(size(bd.Hardware));
%             else
%                 %locMask = ismember(bd.element_info.hardwareType, {'subdural','micro-subdural'});
%                 locMask = ones(size(bd.Hardware));
%             end
%             hw = bd.Hardware(vector(locMask) & vector(hemiMask));
            hw = bd.Hardware;
            bd.plotHardware(hw, varargin{:});
            

        end
    
        function setOpacity(bd, opacity)
            % setOpacity(opacity) sets alpha value of surface to opacity
            % opacity is in [0,1] with 1 being opaque
%             try ax = bd.FigureData.figure.CurrentAxes;
%             catch, ax = gca();
%             end
%             if isempty(ax), error('No axis on braindData figure'); end
            ax = gca;
            patches = findobj(ax, 'Type', 'patch');
            for p = 1 : length(patches)
                patches(p).FaceAlpha = opacity;
            end
            
        end
        
        function plotDepths(bd, varargin)
            % plots all depth electrodes
            %  
            % Optional Input:
            %   same as plot()
            %
            fprintf('Consider using setOpacity function; depths may be hidden within surface\n');
            
            ip = inputParser;
            ip.KeepUnmatched = true;
            ip.addParameter('hemi','', @ischar);
            ip.parse(varargin{:});
            hemi = ip.Results.hemi;
            
            hemiMask = ones(length(bd.Hardware), 1);
            if ~isempty(hemi)
                assert(ismember(hemi, {'lh','rh'}));
                hemiMask = strcmpi(bd.element_info.whichHemi, hemi);
            end
            
            hTypeMask = strcmpi(bd.element_info.hardwareType, 'depth');
            
            hw = bd.Hardware(hTypeMask & hemiMask);
            bd.plotHardware(hw, varargin{:});
        end
        
        function lgd = legend(bd, varargin)
            % Adds legend to current plot. 
            %
            % Optional Input:
            %   axis - handle to axis to make legend for
            %   vertical - pass true to make vertical legend
            %
            % Output:
            %   lgd - legend object. You can change its properties here.
            %
            % Works automatically if you used a
            % hardware-based plot. If you used plotPoint, you need to set the
            % DisplayName property of each point to the appropriate label. See the
            % parameter "legend" of plotPoint
            
            if isempty(bd.FigureData.figure)
                fprintf('Braindata object has no figure handle.\n');
                fprintf('Consider calling plot or (carefully) using bd.FigureData.figure = gcf.\n');
                return;
            end
            
            ip = inputParser;
            ip.addParameter('axis', [] );
            ip.addParameter('vertical', false);
            ip.parse(varargin{:});
            ax = ip.Results.axis;
            isVert = ip.Results.vertical;
            
            default_height = 0.06; % use lgd.Height to change
            
            if isempty(ax)
                fh = bd.FigureData.figure;
                ax = fh.CurrentAxes;
            end
            
            groups = findobj(ax, 'Type', 'hggroup');
            if ~isempty(groups)
                groups = groups(~cellfun('isempty', {groups.DisplayName}));
                groups = groups(~cellfun(@(s) strncmp(s,'data',4), {groups.DisplayName}));
                groups = flip(groups);

                lgd = legend(groups);
                lgd.Interpreter = 'none';
                if ~isVert
                    lgd.Orientation = 'horizontal';
                    lgd.Position = [0, 0, 1, default_height];
                end
            else
                legend(findall(gca, 'type', 'patch'))
            end
            
            
            
        end
        
        function flip(bd) %#ok<MANU>
            % bd.flip() flips brain upside-down 
            % This sets gca.ZDir,ax.Ydir = reverse/normal
            ax = gca;
            switch ax.ZDir
                case 'normal', ax.ZDir = 'reverse'; ax.YDir = 'reverse';
                case 'reverse', ax.ZDir = 'normal'; ax.YDir = 'normal';
            end     
            
        end
        
        function plotDeformSprings(bd, chanNum, fast)
            % plots the edges used in the deformation algorithm
            
            % I assume that leads.csv is in order of increasing chanNum
            % (need to sort if assumption is not valid)
            
            filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, bd.pathEdges);
            assert(exist(filename, 'file') > 0, 'File does not exist');
            load(filename, 'edges');
            edges = sort(edges, 2); %#ok<NODEF> % (low first)
            edges = unique(edges,'rows'); 
            
            subdurals = bd.Jacksheet(strcmp(bd.Jacksheet.hardwareType, 'subdural'), :);
            subdurals = join(subdurals, bd.leads);
            
            if exist('chanNum', 'var') && ~isempty(chanNum)
                if ischar(chanNum)
                    chanNum = find(strcmpi(subdurals.chanName,chanNum));
                end
                edges = edges(edges(:,1)==chanNum | edges(:,2)==chanNum,:);
            else
                % check the edges are on the set of subdurals
                assert(max(edges(:)) == height(subdurals) && min(edges(:)) == 1, 'change code to include hemispherics (maybe others?) in projection');
            end
            
            % plot each line segment
            chan_last = [];
            for i = 1 : length(edges)

                rows = subdurals(edges(i, [1,2]), :);
                chan = rows.chanName{1};
                
                if ~strcmp(chan, chan_last) && (~exist('c','var') || ~exist('chanNum','var'))
                    c = .2 + .8 * rand(1,3);
                end
                
                line(rows.x, rows.y, rows.z, 'color', c, 'lineWidth',2 );
                %bd.plotPoint([rows.x(1), rows.y(1), rows.z(1)], 'color', c);
                %bd.plotPoint([rows.x(2), rows.y(2), rows.z(2)], 'color', c);
                d = norm([rows.x(1), rows.y(1), rows.z(1)] - [rows.x(2), rows.y(2), rows.z(2)]);
                
                chan_last = chan;
                
                
                
                if exist('fast','var') && fast, continue
                else
                    %pause(.1);
                    fprintf('%s --> %s \t %2d cm\n', rows.chanName{1}, rows.chanName{2}, d/10);
                end
            end
        end
        
        function plotAnchors(bd, varargin)
            % Plots anchor points at their actual anchor'd coordinates
            anchors = readtable(fullfile(bd.meta.rootEEGdir, bd.meta.subj, bd.pathAnchors));
            xyz = [anchors.x anchors.y anchors.z];
            bd.plotPoint(xyz, varargin{:});
            
        end
        
        function plotPialProj(bd, chans)
            % Plots projection of electrodes onto standard nodes via tal/intermediates/locs_3_anat/coord2nodes.std141.mat
            filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, 'tal_old/intermediates/locs_3_anat/coord2nodes.std141.mat');
            assert(exist(filename, 'file') > 0, 'Cannot find %s', filename);
            data = load(filename);
            coord2nodes = data.coord2nodes; % chaName, nodes, surfName, hemi
            
            coord2nodes = coord2nodes(strcmpi(coord2nodes.whichHemi, bd.meta.hemi), :);
            fprintf('Using hemi: %s\n', bd.meta.hemi);
            for i=1:length(coord2nodes.nodes)
                bd.plotRegion(coord2nodes.nodes{i}, 'color',[0 1 0]);
            end
            
            
        end
        
        function plotRegion(bd, V, varargin)
            % plotRegion(bd, V, varargin)
            %
            % Colors the set of vertices given by indices in V
            %
            % Optional Parameters:
            %   color - rgb color of region
            %   blend - opacity in [0,1]
            %   boundaryWidth - pass a boundary width. Default no boundary. (2 is good)
            %   useNewSurf - true to create a new surface for the region
            %
            % The bd surface will get colored at all points in V. useNewSurf will
            % instead leave bd.surf unchanged and overlay an translucent new surface.
            % UseBoundary always creates a new surface with just the boundary colored.
            % If you use boundary, appearance is better if you also useNewSurf,
            % although you will be adding 2 new patches which will slow performance if
            % you are plotting a lot.
            %
            
            ip = inputParser;
            ip.addParameter('color',[0 0 1]);
            ip.addParameter('blend',[]);
            ip.addParameter('name','boundary');
            ip.addParameter('boundaryWidth',0);
            ip.addParameter('useNewSurf',0);
            ip.parse(varargin{:});
            color = ip.Results.color;
            blend = ip.Results.blend;
            name = ip.Results.name;
            boundaryWidth = ip.Results.boundaryWidth;
            useNewSurf = ip.Results.useNewSurf;
            
            if ~isempty(bd.FigureData.figure) && ~ishandle(bd.FigureData.figure)
                bd.FigureData.figure = gcf;
            end
            
            
            ax = gca;
            surf = findall(ax, 'type','patch');
            if isempty(surf), error('no surf'); end
            surf = surf(end);
            
            if isempty(blend)
                if useNewSurf,  blend = 1;
                else,           blend = 0.25;
                end
            end
            
            % Add color
            V = unique(V, 'stable');
            if size(color,2) == 3
                % rgb
                fcvd = surf.FaceVertexCData;
                fcvd(V, :) = blend * repmat(color, length(V), 1) + fcvd(V,:) * (1-blend);
                if ~useNewSurf
                    surf.FaceVertexCData = fcvd;
                end
                use_rgb = 1;
            else
                use_rgb = 0;
            end

            
            % Boundary
            if useNewSurf || boundaryWidth > 0
                [Vb,Fb] = bd.getRegionBoundary(V, surf);
                
            
                if boundaryWidth > 0
                    % Only color boundary
                    rgb = nan(length(surf.Vertices), 3);
                    rgb(Vb,:) = repmat(color,length(Vb),1);
                    
                    bnd_surf = patch(...,
                        'faces',Fb,...
                        'vertices',surf.Vertices,...
                        'edgecolor','interp',...
                        'facecolor','interp',...
                        'facevertexcdata',rgb,...
                        'LineWidth', boundaryWidth, ...
                        'specularstrength',0, ...
                        'faceAlpha',blend, ...
                        'faceLighting', 'gouraud', ...
                        'tag','boundary');
                end

                if useNewSurf
                    if use_rgb
                    	CData = repmat(color, length(surf.Vertices), 1);
                    else
                        CData = color;
                    end
                    
                    region_surf = patch(...,
                        'faces',Fb,...
                        'vertices',surf.Vertices,...
                        'edgecolor','none',...
                        'facecolor','interp',...
                        'faceAlpha', blend, ...
                        'facevertexcdata',CData,...
                        'specularstrength',0, ...
                        'faceLighting', 'gouraud', ...
                        'tag','region');
                end
            end
        end
        
        function v = getMedialNodes(bd)
            hem = bd.surf.hemi;
            if isempty(bd.parc) || ~isfield(bd.parc.a2009s,'hemi') || ~strcmpi(bd.parc.a2009s.hemi, hem) || isfield(bd.parc.a2009s,'whichHemi')
                bd.fn_load_parc('parcType','a2009s'); 
            end
            
            
            
            lut = bd.parc.a2009s.lut;
            
            
            
            
            medial_label = sprintf('wm_%s_Unknown',hem);
            v = lut{strcmpi(lut.label_name, medial_label), 'vert_idx'}{1};
            
        end
        
        function sa = getSurfaceArea(bd)
            sa = trimeshSurfaceArea(bd.surf.vertices, bd.surf.faces);
        end
        
        function clearPoints(bd)
            %if ~ishandle(bd.FigureData.figure) bd.FigureData.figure=gcf; end
            
            f = gcf;
            %bd.FigureData.figure
            points = findall(f, 'Tag','electrode');
            for i=1:length(points)
                points(i).Parent = [];
            end
            groups = findall(f, 'Tag', 's');
            for i =1:length(groups)
                groups(i).Parent = [];
            end
            texts = findall(f, 'type', 'text');
            for i = 1:numel(texts)
                texts(i).Parent = [];
            end
        end
        
        function plotRoiData(bd, roi_data, rois)
            
            cm = redWhtBlu(10);
            colormap(cm);
            n = length(cm);
            lims = caxis();
            
            roi_means = nanmean(roi_data,2);
            
            % lay down a base layer
            for i = 1:length(roi_means)
                if ~isnan(roi_means(i))
                    c = [.5 .5 .5];
                    bd.plotRegion(rois{i}, 'color',c, 'blend',1);
                end
            end
            
            
            blend = .5;
            
            
            for i = 1:length(roi_means)
                if ~isnan(roi_means(i))
                    val = roi_means(i);
                    c = cm(1 + round(n * (max(val-lims(1), 0) / (lims(2)-lims(1)))), :);
                    bd.plotRegion(rois{i}, 'color',c, 'blend',blend);
                    %pause(.1)
                end
            end
        end
        
        function vertexSets = roic2roi(bd, roics, rroi_mm)
            % returns sets of mesh indices which fall within rroi_mm of given roics
            %
            % Input:
            %   roics - n-length list of ROICs as *index-in-mesh*
            %   rroi_mm (optional) - geodesic radius to grow ROIs (usually subject-specific)
            % Output:
            %   vertexSets - n-length cell array of mesh indices for each given ROIC
            
            if nargin < 3, rroi_mm = 14; end
            
            
            vertexSets = cell(length(roics), 1);
            
            % load dists if need be
            if isempty(bd.roi) || ~isfield(bd.roi,'dists') || isempty(bd.roi.dists)
                assert(isscalar(bd.meta.hemi) || ischar(bd.meta.hemi), 'Dont know which hemisphere to lead');
                %roi_filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, bd.pathRoi, sprintf('d_roi2std_%s.mat', char(bd.meta.hemi))); 
                roi_filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, 'tal/zloc/roi/roic_Mesh_LUT.csv');
                
                if strfound(bd.surf.dir, 'fsaverage')
                    roi_filename = ['~/ROI/fsaverage/tal/intermediates/roi_1_grow/' sprintf('d_roi2std_%s.mat', char(bd.meta.hemi))];
                end
                %data = load(roi_filename);
                t = readtable(roi_filename);
                
                bd.roi.dists = t(:,1:3); 
                bd.roi.dists.Properties.VariableNames = {'v','d','r'};

                density_mask = ismember(bd.roi.dists.r, bd.roi.roics);
                bd.roi.dists = bd.roi.dists(density_mask,:);
                
                if numel(unique(bd.roi.dists)) < 0.95 * numel(bd.roi.roics)
                    error('There is a suspicously low number of overlapping ROICs between roi_1_grow table and standard ROIC list');
                    % check ~/ROI/fsaverage/tal/intermediates/roi_1_grow/' sprintf('d_roi2std_%s.mat'
                    % against /Users/trottams/eeg_toolbox/localize/other/regions/zroi_base_centers_rh.mat
                end
                
            end
            
            
            % Use roi LUT to find ROI vertex set
            t = bd.roi.dists;   % note v, d, r columns: vertex, distance, center
            d_mask = (t.d <= rroi_mm);
            t = t(d_mask, :);
            for i = 1:length(roics)
                roic_mask = (t.r == roics(i));
                roi_verts = unique(t{roic_mask, 'v'});
                vertexSets{i} = roi_verts;
            end
        end
        
        function [roi_data, findMask] = data2ROIC(bd, data, chanNames, use10mm, rroi_mm)
            % group channel data into rois. 
            %   averages all electrodes within each ROIC's radius
            %
            % Usage: [roi_data, rois, findMask] = data2ROI(bd, data, [chanNames, use10mm, rroi_mm])
            % Outputs
            %   roi_data - matrix of ROIC x channel
            %   findMask - mask into roi_data rows which have any data
            %
            % data should be a vector with 1 value per bd.leads or you can pass specific chanNames
            % roi_data has columns ordered by chanNames, rows ordered by ROIC
            % to get n in region, use sum(~isnan(roi_data),2);
            %
            % example usage:
            %   [rd, mask] = bd.data2ROIC(cgamma, chans, 0, RROI);
            %   roicd = nanmean(rd, 2);
            %   V = bd.roic2roi(bd.roi.roics(mask), RROI);
            %   bd.plotRegionsData(V, roicd(mask), [], [], 1.0, 0);
            
            
            
            USE_MONOPOLAR = 0;
            
            BIPO = 1;
            
            if nargin < 3 || isempty(chanNames), chanNames = bd.leads.chanName; end
            if nargin < 4 || isempty(use10mm), use10mm = 0; end
            if nargin < 5 || isempty(rroi_mm), rroi_mm = 9; end
            
            
            assert(numel(data) >= numel(chanNames));
            
            if BIPO
                temp = bd.meta;
                t = readtable(fullfile(temp.rootEEGdir, temp.subj, 'tal/zloc/mr_pre/coords_mid_euclid.csv'));
                leads = t.chanName;
            else
                leads = bd.leads.chanName;
                assert(isempty(setdiff(chanNames, leads)), 'you gave chanNames not in leads');
            end
            
            if USE_MONOPOLAR
                assert(all(cellfun(@strcmp,bd.leads.chanName,{bd.Monopolar.chanName}')));
            end
            
            % given chanNames index into leads; thus i'th data point must correspond to i'th chanNdx
            chanNdx = cellfun(@(chan) find(strcmpi(chan,leads),1), chanNames);
            
            bd.setRoics(use10mm);
            
            % for every roic we may have multiple data points
            roi_data = nan(numel(bd.roi.roics), numel(data));
            %rois = cell(numel(bd.roi.roics), 1);
            
            if ~USE_MONOPOLAR
                if BIPO
                    xyz = t{:,{'x','y','z'}};
                    loaded = load(fullfile(temp.rootEEGdir, temp.subj, 'tal/zloc/roi/ROIC_lead_LUT_bipolar.mat'));
                    t = loaded.(char(fieldnames(loaded)));
                else
                    xyz = bd.leads{:,{'x','y','z'}};
                    loaded = readtable(fullfile(temp.rootEEGdir, temp.subj, 'tal/zloc/roi/ROIC_lead_LUT.mat'));
                end
                
                %custom_std_nodes = lead_to_pial(xyz, [], bd.surf.vertices, 1.5, 0);
                custom_std_nodes = t.lead_verts;
            end
              
            for i = 1:length(chanNames)
                % 1) get nearest pial mesh notes
                if USE_MONOPOLAR
                    mesh_nodes = bd.Monopolar(chanNdx(i)).std_nodes;
                else
                    mesh_nodes = custom_std_nodes{chanNdx(i)};
                end
                [roics, ~, roicNdx] = bd.vertex2ROI(mesh_nodes, use10mm, rroi_mm);
                roi_data(roicNdx,i) = data(i);
                %rois(roicNdx) = repmat({roi},length(roicNdx),1);
                %bd.plotPoint(std_nodes(1),'color',[0 0 0]);
            end
            % The above loop results in a sparse matrix where data is put into columns of each ROIC row
            
            % mask of any ROICs with datacount > 0;
            findMask = logical(sum(~isnan(roi_data),2) > 0);
            
            
        end
        
        function [roics, roi, roicNdx] = vertex2ROI(bd, vertex, use10mm, rroi_mm, roi_filename)
            % Input:
            %   vertex      - vertex set to map
            %
            % Input (Optional):
            %   use10mm     - 1 [default] for 10mm inter-ROIC spacing, 0 for 5mm
            %   rroi_mm     - ROI radius (default Infinite)
            %   roi_filename- 
            % Output:
            %   roics       - roi centers as surface vertices
            %   roi         - vertices within rroi_mm of any ROI center
            %   roicNdx     - index into bd.roi.roics of roi centers
            
            USE_MONOPOLAR = 0;
            
            if nargin < 3 || isempty(use10mm), use10mm = 0; end
            if nargin < 4 || isempty(rroi_mm), rroi_mm = Inf; end
            if nargin < 5 || isempty(roi_filename), roi_filename = []; end
            
            
            %bd.roicDists=[];   
            
            if isempty(roi_filename) && USE_MONOPOLAR
                fprintf('Using monopolar.mat\n');
                % use monopolar.mat
                mp = bd.Monopolar;
                
                switch use10mm
                    case 1
                        rois = mp.zroi_std_10mm;
                        dist = mp.d_elec2_10mm;
                    case 0
                        rois = mp.zroi_std_5mm;
                        dist = mp.d_elec2_5mm;
                end
            
                d_mask = dist <= rroi_mm;
                
                
            else
                %fprintf('Loading from file instead of using monopolar.mat\n');
                % Load from 'tal/intermediates/roi_1_grow/'
                
                bd.setRoics(use10mm);
                
                % Get subject specific ROI growth for roics
                if isempty(bd.roi) || ~isfield(bd.roi,'dists') || isempty(bd.roi.dists)
                    assert(isscalar(bd.meta.hemi) || ischar(bd.meta.hemi), 'Dont know which hemisphere to load');
                    if isempty(roi_filename)
                        roi_filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, bd.pathRoi, sprintf('d_roi2std_%s.mat', char(bd.meta.hemi))); 
                    end
                    if exist(roi_filename,'file')
                        data = load(roi_filename);
                        bd.roi.dists = data.d_roi2std(:,1:3); 
                    else
                        roi_filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, 'tal/zloc/roi/ROIC_Mesh_LUT.csv');
                        t = readtable(roi_filename);
                        bd.roi.dists = t(:,{'vertex','d','ROIC_mesh_ndx'});
                    end
                    
                    %% 
                    bd.roi.dists.Properties.VariableNames = {'v','d','r'};

                    density_mask = ismember(bd.roi.dists.r, bd.roi.roics);
                    bd.roi.dists = bd.roi.dists(density_mask,:);
                    %% 
                end
                
                
                %                               (1)      (2)
                % general procedure is electrode --> ROIC --> ROI vertices
                t = bd.roi.dists;   % note v, d, r columns: vertex, distance, center
                d_mask = t.d <= rroi_mm;
                v_mask = ismember(t.v, vertex);  
                roics = unique(t.r(d_mask & v_mask));   % 1 roics which are within distance of query vertices
                c_mask = ismember(t.r, roics);
                roi = unique(t.v(d_mask & c_mask));       % 2 vertices which are within distance of centers
                roicNdx = arrayfun(@(x) find(x==bd.roi.roics), roics,'uniformOutput',0);
                roicNdx = cell2mat(roicNdx(~cellfun('isempty',roicNdx)));
            end
            
            
        end

        function all_roics = show_elect2roi(bd, rroi_mm, varargin)
            % This will store each pair in two columns: 
            %   1 - index into monopolar
            %   2 - index into mesh;

            %monopolar = bd.Monopolar;
            %monopolar = monopolar(strcmpi({monopolar.loc_hemi}, bd.meta.hemi));
            
            roi_filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, 'tal/zloc/roi/ROIC_lead_LUT.mat');
            data = load(roi_filename);
            lead_roi_lut = data.lead_roi_lut;
            
            all_roics = [];
            
            WAIT = 1;
            
            % load roi dist table
            bd.roic2roi([], rroi_mm);
            
            % From nearest pial nodes, look up centers within growth distance
            nodes_per_e = lead_roi_lut.ROIC_verts;
            mind = 0;
            for i = 1 : length(nodes_per_e)
                enodes = nodes_per_e{i};

                % filter LUT to rows with our electrode's nodes
                rows = bd.roi.dists(ismember(bd.roi.dists.v, enodes), :);
                mind = min(rows.d);

                % Find ROICs that contain this electrode in their r_growth ROI
                mask = rows.d < rroi_mm;
                rows = rows(mask, :);
                roics = unique(rows.r);

                if isempty(roics)
                    fprintf('%s %s %s is not in an ROI...', bd.meta.subj, bd.meta.hem, monopolar(i).chanName);
                    fprintf('closest roi is %.03f mm away\n', mind);
                else
                    bd.plotPoint(roics, 'radius',0.75);
                    all_roics = union(all_roics, roics);
                end
                
                switch WAIT
                    case 1
                        pause(.1);
                    case 2
                        keyboard;
                end
            end
        end
        
        function lights = camlight(bd,intensity)
            % many camlights
            if nargin < 2, intensity = 2.8; end
            lights = findall(gca,'type','light');
            for i=1:length(lights), lights(i).Parent=[]; end
            lights = {};
            lights = [lights {camlight(0,0)}];
            lights = [lights {camlight(90,0)}];
            lights = [lights {camlight(0,90)}];
            lights = [lights {camlight(90,90)}];
            lights = [lights {camlight(-90,0)}];
            lights = [lights {camlight(0,-90)}];
            lights = [lights {camlight(-90,-90)}];
            lights = [lights {camlight(-90,90)}];
            for i=1:length(lights)
                lights{i}.Color = [1 1 1] ./ length(lights) .* intensity;
                %lights{i}.Style = 'infinite';
            end
        end
        
        function clear(bd, varargin)
            %f = bd.FigureData.figure;
            %if ~ishandle(f), f = gcf; end
            f = gcf;
            clf(f);
            bd.plot('handle',f, varargin{:});
        end
        
        function [vertexData, vertexCnt] = plotRegionsData(bd, vertexSets, data, clim, cmap, blend, cumulative)
            % plots data on multiple vertex sets
            % Input:
            %   vertexSets  - n x 1 cell array of surface index arrays
            %   data        - n x 1 array of data values to plot
            %
            % Optional Input:
            %   clim        - [low high]
            %   cmap        - colormap
            %   blend       - [default=1] 1 gives full color, 0.1 gives 10 percent color
            %   cumulative  - if 0 [default], vertex "x" is divided by m=(# datapoints at x)
            %                 if 1, vertex "x" is sum_{datapoints at x}(data) at x
            % Output:
            %   vertexData  - data plotted at each vertex from [1:max(vertexSets)]
            %   vertexCnt   - how many data points are at each vertex
            
            if nargin < 4 || isempty(clim), clim = caxis(); end
            if nargin < 5 || isempty(cmap), cmap = colormap(); end
            if nargin < 6 || isempty(blend), blend = 1; end
            if nargin < 7 || isempty(cumulative), cumulative = 0; end
            
            n = length(vertexSets);
            assert(length(data) == n, 'data and vertexSets must be the same size');
            
            
            hiVertex = cellfun(@max, vertexSets);
            hiVertex = max(hiVertex);
            
            % 1 row for every vertex (sparse)
            vertexData = zeros(bd.stdNumVert, 1);
            vertexCnt  = zeros(bd.stdNumVert, 1);
            
            
            % aggregate data from different vertex sets onto each vertex
            for i = 1:n
                vertexCnt(vertexSets{i}) = vertexCnt(vertexSets{i}) + 1;
                vertexData(vertexSets{i}) = vertexData(vertexSets{i}) + data(i);
            end
            
            % average at each vertex
            findMask = vertexCnt > 0;
            vertexData(~findMask) = nan;
            if ~cumulative
                vertexData(findMask) = vertexData(findMask) ./ vertexCnt(findMask);
            end
            
            

            
            bd.plotRegion(find(findMask),'usenewsurf',1,'color',vertexData, 'blend',blend);
            return
            
            % map onto colormap
            val2col = @(val) cmap(min(max(round(1+ (length(cmap)-1) * ( (val-clim(1)) ./ (clim(2)-clim(1)) )), 1),length(cmap)), :);
            vertexColor = nan(sum(findMask), 3);
            vertexColor = arrayfun(val2col, vertexData(findMask), 'uniformOutput',0);
            
            % get surf
            try surf = findall(gcf,'type','patch'); surf=surf(1);
            catch e, surf=[]; end
            
            % color
            %rgb = surf.FaceVertexCData;
            %rgb(findMask,:) = repmat(base,sum(findMask),1);
            surf.FaceVertexCData(findMask,:) = blend * cell2mat(vertexColor) + (1-blend) * surf.FaceVertexCData(findMask,:);
            
            
            
        end
        
        function setRoics(bd, use10mm)
            if nargin < 2, use10mm = 0; end
            
            % set roics to be a subset of all the centers (eg 600 out of 2400)
            if isempty(bd.roi) || ~isfield(bd.roi, 'roics') || isempty(bd.roi.roics) || (use10mm && bd.roi.density==5) || (~use10mm && bd.roi.density==10)
                one_hemi = bd.meta.hemi;
                if iscellstr(one_hemi), one_hemi = one_hemi{1}; end
                roic_filename = which(sprintf('zroi_base_centers_%s.mat',char(one_hemi)));
                %fprintf('ROIC file: %s\n', roic_filename);
                data = load(roic_filename);
                all_centers = data.vert_idx;
                % 2400 = 5mm, 600=10mm
                if use10mm
                    bd.roi.density = 10;
                    bd.roi.roics = all_centers(1:bd.ncenter_course);
                else
                    bd.roi.density = 5;
                    bd.roi.roics = all_centers(1:bd.ncenter_fine);
                end

            end
            
            
        end
        
        function t = getXyzBipolar(bd)
            % loads the biplar xyz table (euclidean) which has columns
            %   chanName (e.g. G1-G2)
            %   x, y, z
            %   whichHemi
            %   hardwareType
            fname_temp_bipolar = fullfile(bd.meta.rootEEGdir, bd.meta.subj, 'tal/zloc/mr_pre/coords_mid_euclid.csv');
            if exist(fname_temp_bipolar,'file')
                t = readtableSafe(fname_temp_bipolar);
            end
            
        end
        
    end % methdos public
    
    methods (Access = private)
        function nodes = getChanNodes(bd, chanNames, single)
        % chanName --> vertices
            if nargin > 2 && single
                firstNode = true;
            else
                firstNode = false;
            end
        
            % intermediate
            if isempty(bd.Monopolar)
                filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, 'tal/zloc/roi/ROIC_lead_LUT_monopolar.mat');
                %filename = fullfile(bd.meta.rootEEGdir, bd.meta.subj, 'tal/zloc/roi/ROIC_lead_LUT_bipolar.mat');
                data = load(filename); fn = fieldnames(data);
                coord2nodes = data.(fn{1});
                nodes = coord2nodes(ismember(coord2nodes.chanName,chanNames), :).lead_verts; %#ok<NODEF>
            else
                % monopolar
                chans = {bd.Monopolar.chanName};
                mask = ismember(chans, chanNames);
                nodes = {bd.Monopolar(mask).std_nodes};     
            end
            
            % non empty
            nodes = nodes(~cellfun('isempty',nodes));

            % nothing returned
            if isempty(nodes)
                fprintf('Error getting chanNodes: no nodes found for channels\n');
            end
            
            % get only first node
            if firstNode
                nodes = cellfun(@(x) x(1), nodes);
            end
            
            % dereference scalar
            %if length(nodes) == 1, nodes = nodes{1}; end
        end
        
        function surf = getASurf(bd)
        % gets one of the lh/rh surfaces. undefined which one it gives you
            if isfield(bd.surf, 'lh'), surf = bd.surf.lh;
            elseif isfield(bd.surf, 'rh'), surf = bd.surf.rh;
            else surf = bd.surf;
            end
        end
        
        function [V,F] = getRegionBoundary(bd, V, surf)
            % Gets bounding edges, returned in order
            F = surf.Faces;
            
            % Get faces of region nodes
            V2F = surfing_nodeidxs2faceidxs(F');
            Fndx = V2F(V,:);
            Fndx = unique(Fndx(:));
            Fndx = Fndx(Fndx~=0); % V2F had a zero row and column
            %Fr = unique(sort(F(Fndx,:), 2), 'rows');
            F = unique(F(Fndx,:), 'rows', 'stable');
            % Fr is now the faces of the region's nodes
            
            % Edges
            E = [F(:,1) F(:,2); F(:,2) F(:,3); F(:,1) F(:,3)];
            E = sort(E, 2);
            [Eu,ie] = unique(E, 'rows', 'stable');
            i_not_unique = setdiff(1:length(E), ie);
            E_mult = E(i_not_unique, :);            % E with 1 copy of every edge remvoed
            E_bnd = setdiff(Eu, E_mult, 'rows');    % Unique edges without edges that appear in >1 face

            V = unique(E_bnd(:));
        end
        
        
        
    end
    
    methods (Static)
        % Apparently there is no freesurfer/std.141.a2009s.annot file
%         function v = getFSAverageMedial(hem)
%             bd = braindata;
%             bd.loadFSAverage('hemi',hem);
%             bd.meta.suma_dir = '/Users/trottams/Applications/freesurfer/subjects/fsaverage/SUMA';
%             bd.fn_load_parc('hemi',hem,'parcType','a2009s');
%             lut = bd.parc.a2009s.lut;
%             medial_label = sprintf('wm_%s_Unknown',hem);
%             v = lut{strcmpi(lut.label_name, medial_label), 'vert_idx'}{1};
%             
%         end
        
        

    end
end % class

%         function parcellate_lobes(bd)
%             bd.fn_load_parc(
%             
%             function fn_load_parc(bd, varargin)
%         % FN_LOAD_PARC does ____
%         % must call with surface loaded as combined
%         %
%         % Accepts the following (optional) key-value pairs:
%         %   subj - will use class property if not passed
%         %   hemi - uses surface's property if not overriden
%         %   surfRez - uses surface's property if not overriden
%         %   parcType [a2009s]
%         %   dir
%             try
%                 default_cell = {
%                     'subj',bd.meta.subj,'parameter';
%                     'hemi', bd.surf.hemi, 'parameter';
%                     'surfRez',bd.surf.surfRez,'parameter';
%                     'parcType','a2009s','parameter';
%                     'util2dir',bd.meta.util2dir,'parameter';
%                     'dir',bd.meta.suma_dir,'parameter'
%                     };
%             catch e
%                 warning('You may have tried to call load_parc without using a combined surface...');
%                 fprintf('Caught this error: %s\n', e.message);
%                 fprintf('I will relad the surface with combined and try it again for you.\n');
%                 keyboard;
%                 fprintf('Calling fn_load_surface("hemi","combined..."\n');
%                 bd.fn_load_surf('hemi','combined');
%                 bd.fn_load_parc(varargin{:});
%                 return;
%             end
%             inputs = util_extract_inputs(default_cell, varargin);
%             
%             % now load the parc file
%             bd.parc.(inputs.parcType) = surf_load_parc(inputs);
%         end
%         end