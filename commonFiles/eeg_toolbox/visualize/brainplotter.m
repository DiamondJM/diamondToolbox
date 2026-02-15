%   BrainPlotter is responsible for:
%       - Storing surface geometries
%       - Data sets (underlays like curvature, sulc, atlases, + overlays like user data)
%       - XYZ electrode points
%   
%   Brainplotter is meant to accompany braindata2. Previously, plotting functions and subject-specific data-loading
%       were all done in one braindata object. Now braindata2 just does data-loading and storing, and plotting and figure
%       functions are done in brainplotter.
%
%   brainplotter Properties:
%       surfaces        - A struct which stores a set of all loaded surfaces
%
%   brainplotter Methods:
%       loadSurface     - Loads surface by directory or filepath
%       unload          - Unloads all surfaces
%       chooseNewSurface - Provides a menu of all possible surfaces you can load
%
%       plot            - Renders a given surface on the current axis
%       plotPoint       - Plots a given xyz-coordinate (or multiple points)
%       plotPointGroup  - Plots a set of points and color codes them by group
%       unplot          - Removes brain surfaces from current axis
%       clearPoints     - Removes all plotted points from current axis
%       
%       plotRegion      - Shade a given region with a given color
%       plotRegionsData - Shade a region according to a dataset
%       plotSulc        - Shade the surface by gyrus/sulcus binary category
%
%       plotAtlas       - Colors the surface according to an atlas file (Note: we only have SUMA data for destrieux atlas currently)
%
%       viewAll         - Utility for saving figures; copyies current axis' objects to new figure with standard rotations
%       saveToRGB       - RGB-ifies a given Indexed-colored matlab patch. (See help(plotRegion) for details on patch color modes)
%       
%
%   More help:
%       - Create an instance of brainplotter by calling "bp = brainplotter()".
%         (if you have a braindata2 object called "bd", you'll likely want to next call "bp.loadSurface(bd)" )
%       - You can use the MATLAB functions "methods(bp)" and "properties(bp)" to investigate your object
%       - Try typing: "doc brainplotter"
%       - See the tutorial on /Volumes/Shares/FRNU/labmeetings/2018_03_20_trotta_braindata/bdtutorial/braindata_tutorial.m
%           It has a FAQ and a walk-through. Please copy it locally before running it.
%   
% See also: braindata2, braindata

% REVISION HISTORY 
%   ( ** Note: remember to update the version() function! ** )
%   04/2018 MST     - Use getToolboxDir instead of hard-coded path
%                   - Add replotLabels() function
%   05/2018 MST     - Fix bug where plotting points with RGB color changed the colorbar CLim
%                   - Updated camlights positions to be more uniform + surface color tweaks
%                   - saveToRGB now public and preserves nans
%   08/2018 MST     - Fix incorrect label bug when plotting with 'useNode',1
%                   - Add plotAtlas (Note: we only have SUMA data for destrieux atlas currently)

classdef brainplotter < handle
    properties
        surfaces = struct(); % Struct of surfaces
    end
    
    properties (Constant) % theoretically private, public for development convenience
        stdNumVert  = 198812; % per hemisphere
        dataPath    = fullfile(getToolboxDir(), 'visualize/data');
    end
    
    properties %(private) public for development convenience
        rngseed
    end
    
    methods % Constructor (and other standard class functions)
        function bp = brainplotter()
            if ~exist(brainplotter.dataPath, 'dir')
                warning('You need to set the constant dataPath property. File not found: %s', brainplotter.dataPath);
            end
            bp.rngseed = 1;
        end
                
        function disp(bp)
            % Displays the brainplotter in an informative way
            bp.line(0); bp.line(1); bp.line(0);
            
            fprintf('Help:\n\tBrainplotter is a class for loading, storing, and displaying brain surfaces and electrode data.\n\tTo use, call bp = brainplotter and see help(brainplotter).\n\n');
            
            % display all loaded surfaces
            fprintf('Data:\n');
            fns = fieldnames(bp.surfaces);
            fprintf('\tLoaded bp.surfaces fields:\n');
            if numel(fns) == 0
                fprintf('\tNo surfaces loaded. Try "help bp.loadSurface"\n');
            end
            for i = 1:numel(fns)
                fprintf('\t\t%s', fns{i});
            end
            fprintf('\n\n');
            bp.version();
            bp.line(0);
        end
        
    end
    
    methods % Loading (and unloading)
        
        function [gsurfs,names] = loadSurface(bp, fname, name)
            % loadSurface(bp, fname, name) Loads a gifti surface by directory or filepath
            %
            % Inputs
            %   fname   - name of directory or full filepath to gifti surface (usually SUMA/). 
            %             If not passed, loads the average brain.
            %             If you pass a subject-loaded BRAINDATA2 object, will use that subject's SUMA/ directory by default.
            %             Can also be a cell string of fnames
            %
            %   name    - If fname is a filename, 
            %               loaded surface will be stored in bp.surfaces.(name). Suggested name format: "pial_rh"
            %               If fname was a cell string, pass an equal number of names.
            %             If fname is a directory or BRAINDATA2 object, pass a space or _ delimited list of keywords to filter
            %               the list of available surfaces, (e.g. "rh_pial" loads surfaces which contain "rh" and "pial" in their name)
            %
            % OUTPUTS:
            %   gsurfs
            %   surfNames
            %
            %  Examples:
            %   1) bd=braindata2(subj,root); bp.loadSurface(bd, 'lh pial');
            %   2) bp.loadSurface(pathToGiftiSurfaceFile, 'mySurfaceName')
            
            gsurfs = [];
            names = [];
            
            if nargin < 2 || isempty(fname)
                fprintf(['You gave no file/directory to loadSurface.' ...
			 ' Loading fsaverage by default\n']);

                fname = fullfile(bp.dataPath, 'fsaverage/SUMA');
            end
            if nargin < 3
                name = [];
            end

            
            % Load each fname separately
            if iscellstr(fname)
                k = numel(fname);
                if iscellstr(name) && numel(name) == k
                    names = name;
                else
                    names = [];
                end
                for i = 1:k
                    if ~isempty(names)
                        [gsurfs{i},names{i}] = bp.loadSurface(fname{i}, names{i});
                    else
                        [gsurfs{i},names{i}] = bp.loadSurface(fname{i});
                    end
                end
                return
            end
            
            % braindata2 object passed
            if isa(fname, 'braindata2')
                bd = fname;
                fname = bd.filepaths.surfaces;
                
                % set rngseed based on subject number
                [~,subjNum] = util_split_stringnum(bd.subj);
                if ~isempty(subjNum)
                    bp.rngseed = subjNum;
                end
                
                if exist(fname, 'file')
                    [gsurfs,names] = bp.loadSurface(fname, name);
                    return;
                else
                    fprintf('braindata2.loadSurface: Passed braindata2 object but could not find %s\n', fname);
                end
            elseif isdir(fname)
                [gsurfs,names] = bp.chooseNewSurface(fname, name);
                return;
            elseif ~exist(fname, 'file')
                
                fprintf('Cannot find name given, looking in present directory... (name given: %s)\n', fname);
                
                % maybe they gave a name like lh pial and are looking in the pwd
                [gsurfs,names] = bp.chooseNewSurface(pwd, fname);
                if ~isempty(gsurfs)
                    return
                end
                    
            end
            
            surf_index = 1 + numel(fieldnames(bp.surfaces));
            if nargin < 3 || isempty(name)
                name = sprintf('surface_%d', surf_index); 
            end
            
            assert(exist(fname, 'file') > 0, 'file not found: %s\n', fname);
            
            fprintf('Loading [%s] into [bp.surfaces.%s]...', fname, name);
            gsurfs = gifti(fname);
            names = name;
            bp.surfaces.(name) = gsurfs;
            fprintf('done!\n');
        end
        
        function unload(bp)
            % Unloads all surfaces
            bp.surfaces = struct();
        end
        
        function [gsurfs, names] = chooseNewSurface(bp, d, name)
            % chooseNewSurface(bp, d) Presents a menu to the user to select a surface and calls loadSurface
            %
            % Input:
            %   d       - directory of surfaces
            %   name    - a "_" or space delimited list of keywords to filter surfaces (e.g. "rh pial")
            %
            % OUTPUT:
            %   gsurfs  - cell array of gifti surface structs with properties:
            %               faces
            %               vertices
            %   names   - cell array of surface names
            %
            % Options are all std.141 gifti's in the given directory (or the fsaverage directory by default)
            
            gsurfs = [];
            names = [];
            
            if nargin < 2 || isempty(d)
                fprintf('You gave no file/directory to chooseNewSurface. Loading fsaverage by default\n');
                d = fullfile(bp.dataPath, 'fsaverage/SUMA');
            end
            if nargin < 3
                name = [];
            end
            
            assert(exist(d, 'dir') > 0, 'Directory not found: %s\n', d);
            
            contents = lsCell(d);
            contents = contents(strfound(contents, '.gii') & strfound(contents, 'std.141'));
            if ~isempty(name)
                % Split names
                if contains(name, ' ')
                    name_split = strsplit(name, ' ');
                elseif contains(name, '_')
                    name_split = strsplit(name, '_');
                elseif contains(name, '.')
                    name_split = strsplit(name, '.');
                else
                    name_split = cellstr(name);
                end
                for i = 1:numel(name_split)
                    contents = contents(strfound(lower(contents), lower(name_split{i})));
                end
            end
            
            if isempty(contents)
                error('Nothng found by name %s (also checked in directory %s)', name, d);
            end
            
            if numel(contents) > 1
                bp.hint('Skip the menu by calling loadSurface with an argument that uniquely identifies a surface (like "pial_lh" or "lh pial")');
                title    = sprintf('Choose a surface from %s', d);
                ndx      = menuText(title, contents{:}, 'multiSelect',1);
            else
                ndx      = 1;
            end
            
            for i=1:numel(ndx)
                fname    = fullfile(d, contents{ndx(i)});
                parts    = strsplit(contents{ndx(i)}, '.');
                hemi     = parts{end-2};
                name     = sprintf('%s_%s', parts{end-1}, hemi);

                [gsurf,name] = bp.loadSurface(fname, name);
                gsurfs = cat(1, gsurfs, gsurf);
                names = cat(1, names, name);
            end
        end
        
        function [surfs,names] = get_gsurfs(bp, names)
            % surfs = get_gsurfs(bp, name) 
            % Returns a cell array of loaded gifti surfaces based on name
            
            if nargin < 2 || isempty(names)
                if numel(fieldnames(bp.surfaces)) == 1
                    names = char(fieldnames(bp.surfaces));
                elseif numel(fieldnames(bp.surfaces)) < 1
                    error('No surfaces loaded');
                else
                    [~,names] = menuText('Choose surface(s)',fieldnames(bp.surfaces));
                end
            end
            
            if ~iscellstr(names)
                names = cellstr(names);
            end
            
            % Return each named, loaded surface (if it really has been loaded)
            fns = fieldnames(bp.surfaces); % field names
            surfs = [];
            for i = 1 : numel(names)
                name = names{i};
                if ismember(name, fns)
                    surfs{i} = bp.surfaces.(name);
                else
                    candidate_surfs = fns(contains(lower(fns), lower(name)));
                    if numel(candidate_surfs) == 1
                        surfs{i} = bp.surfaces.(candidate_surfs{1});
                    elseif numel(candidate_surfs) == 0
                        error('No surf loaded matching name %s', name);
                    else
                        % present menu to resolve ambiguity
                        [~,surf_name] = menuText('Choose surface',candidate_surfs, 0);
                        surfs{i} = bp.surfaces.(surf_name);
                    end
                end
            end
            
        end
    end
    
    methods % Plotting
        function [psurf, surfname] = plot(bp, surfname)
            % psurf = plotSurface(bp, name)
            % Plots a loaded surface as a patch. (You must have called loadSurface to load a gifti surface from file into memory)
            %
            % INPUT
            %   surfname - Name of loaded surface. Optional if only 1 surface is loaded
            %   
            % OUTPUT
            %   psurf    - handle to rendered matlab patch object
            %   surfname - name of surface that was resolved/passed-in 
            %
            % DETAILS
            %   After the surface you intend to plot is resolved, the following core functions are called:
            %     1) viz_surf_main.m is called. This is an ancient in-house function which implements color-mixing
            %        and eventually calls matlab's patch function
            %     2) bp.apply_axis_settings() - standard brainplotting axis settings
            %     3) bp.apply_surf_settings() - standard brainplotting patch setttings (e.g. camlights)
            %
            % See Also: loadSurface

            if nargin < 2, surfname = ''; end
            if iscellstr(surfname)
                if numel(surfname) == 1
                    surfname = char(surfname);
                elseif numel(surfname) > 1
                    for i=1:numel(surfname)
                        [psurf{i}, surfname{i}] = plot(bp, surfname{i});
                    end
                    return;
                end
            end
            
            if isempty(surfname) && numel(fieldnames(bp.surfaces)) > 1
                bp.hint('Skip the menu! If you have more than 1 surface loaded, pass the surfname (or cell of surfnames) to the plot function');
            end
            
            % Get gifti surface and plot it
            [surfs, surfname] = bp.get_gsurfs(surfname);
            for i = 1:numel(surfs)
                bp.unplot(surfname{i});
                psurf = viz_surf_main(surfs{i}, 'tag', surfname{i});
            end

            % standard axis settings
            bp.apply_axis_settings();
            bp.apply_surf_settings();
        end
        function clear(bp)
            % clear 
            
            
            regions = findall(gca, 'Tag', 'region');
            for i=1:numel(regions)
                regions(i).Parent = [];
            end
            
            % Reset color to gray
            fn = cellstr(fieldnames(bp.surfaces));
            for i=1:numel(fn)
                name = fn{i};
                psurf = bp.get_psurf(name);
                if ~isempty(psurf)
                    set(psurf, 'FaceVertexCData', repmat(0.5, length(psurf.Vertices), 3))
                end
            end
            
            % remove all other surfaces
            to_remove = findall(gca, 'type', 'patch');
            to_remove = to_remove(arrayfun(@(x) isprop(x, 'Tag'), to_remove));
            if numel(to_remove) > 0
                arrayfun(@(x) set(x,'Parent',[]), to_remove(~ismember({to_remove.Tag},fn)));
            end
            
            % remove text
            texts = findall(gca, 'type', 'Text');
            for i = 1:numel(texts)
                texts(i).Parent = [];
            end
            
            % remove poitns
            brainplotter.clearPoints();
            
            % remove legend
            l = legend();
            delete(l);
        end
        function clearRegions(bp)
            % Clears regions and boundaries
            
            % first clear all non-surface patches
            patches = findall(gca, 'type', 'patch');
            if isempty(patches); return; end
            tags = {patches.Tag};
            pdelete = patches(~ismember(tags, fieldnames(bp.surfaces)));
            puncolor= patches(ismember(tags, fieldnames(bp.surfaces)));
            for i=1:numel(pdelete)
                set(pdelete(i),'Parent',[]);
            end
            
            % Now restore surfaces to their gray color
            for i=1:numel(puncolor)
                set(puncolor(i), 'FaceVertexCData', repmat(0.5, bp.stdNumVert, 3));
            end
            
        end
        function unplot(bp, name)
            % unplot(bp, name) Unplots a surface named "name"
            if nargin < 2, name = []; end
            if ischar(name) name = cellstr(name); end
            
            for i = 1:numel(name)
                surf = bp.get_psurf(name{i});
                if ~isempty(surf)
                    surf.Parent = [];
                end
            end
        end
        function [hgroup, hsphere] = plotPoint(bp, points, varargin)
        % Plots points on the current figure
        %
        % If plotting a subject's point on their own surface, you can pass XYZ
        % coordinates. If plotting onto an avergae or other standard surface, pass
        % indices.
        %
        %
        %
        % Inputs:
        %   points - This can be a variety of things:
        %            - A matrix of xyz coords
        %            - An array of mesh indices
        %            - A cell array of mesh index lists
        %            - A table meeting one of these requirements:
        %               - x, y, and z columns OR
        %               - nearest_mesh_ind column (see note 1)
        %                 * If both x,y,z present and nearest_mesh_ind, pass in
        %                   the UseNode parameter
        %
        % Optional Input
        %   surfname - Name of loaded surface
        %
        % Optional (key-value) Inputs:
        %   plot3   - if 1, use plot3 instead of spherical patch. Default 0
        %   color   - default [0 1 0] green for each point
        %   radius  - default 1.5 (mm)
        %   legend  - string label that applies to all the points passed.
        %               Sets the display name of the group 
        %   ctNum   - If you pass in 1 or 2, when looking up xyz-coords for channel
        %               names, will look in CT_[ctNum]/leads.csv instead of
        %               tal/leads.csv
        %   label   - Plot the text labels given in a cell array
        %   labelColor - Color of text labels
        %   alpha   - alpha level
        %   data    - mapped into colormap
        %   useNode - If there is an option, useNode instead of xyz (default false)
        %   lines   - passing the expected inter-electrode distance (mm) will draw lines at lines +/- 20% inter-elec (see code to adjust)
        %   element_info - if you want to draw the lines with the actual neighbors, pass subject's element info file path or table here
        %   tagName - if (1) you are drawing lines, (2) points is not a table with same-tag chanNames or tagName, and (3) you passed element_info, then 
        %               you should also pass this. If you pass this without points being a table, your points had better be in numerical order
        %   shape   - TODO 'sphere' (default) or 'cylinder'
        %   mesh_clust_d - If you are passing mesh index lists, by default a point is placed at the first mesh index in the list. However, for sulcus-
        %             spanning electrodes, we may want to plot 2+ points at each group of mesh indices. Pass a distance here if you want to clusterize
        %             such that mesh_indices >=~ mesh_clust_d apart will be plotted as separate points. Default is 3. 0 to just use the first index.
        %
        % Outputs:
        %   hgroup  - Handle to group of the point sphere(s)
        %   hsphere - Handle(s) to point sphere(s)
        %
        % Notes:
        %   (1) If you pass a table with nearest_mesh_ind column, this column's cells should be lists of mesh indices. A point will be placed to 
        %       represent each cluster of mesh locations
        %   
        %
        %
        %
        
            mesh_ind = [];
            
            %% ------ argument handling ------ %%
            if istable(points)
                tpoints = points;
                n = height(tpoints);
            else
                n = size(points,1);
                tpoints = [];
            end

            if n==0, return; end

            ip = inputParser;
            ip.addParameter('plot3', 0);
            ip.addParameter('color', repmat([0 1 0], n, 1));
            ip.addParameter('radius', 1.5, @isnumeric);
            ip.addParameter('legend', []);
            ip.addParameter('ctNum', [], @isnumeric);
            ip.addParameter('label', '');
            ip.addParameter('labelColor', [0 0 0]);
            ip.addParameter('alpha',1);
            ip.addParameter('data',[]);
            ip.addParameter('lines',0);
            ip.addParameter('useNode',0);
            ip.addParameter('element_info',[]);
            ip.addParameter('tagName',[]);
            ip.addParameter('shape','sphere');
            ip.addParameter('mesh_clust_d', 3);
            ip.addOptional('surfName', '', @(s) ~ismember(s, ip.Parameters));
            ip.KeepUnmatched = true;
            ip.parse(varargin{:});
            surfName = ip.Results.surfName;
            usePlot3 = ip.Results.plot3;
            color   = ip.Results.color;
            radius  = ip.Results.radius;
            legendLabel = char(ip.Results.legend);
            ctNum   = ip.Results.ctNum;
            label   = ip.Results.label;
            labelColor= ip.Results.labelColor;
            alpha   = ip.Results.alpha;
            data    = ip.Results.data;
            lines   = ip.Results.lines;
            element_info = ip.Results.element_info;
            tagName = ip.Results.tagName;
            shape  = ip.Results.shape;
            mesh_clust_d = ip.Results.mesh_clust_d;
            useNode = ip.Results.useNode;
            
            
            % Type of points
            if ~isempty(tpoints)
                if ismember('x', tpoints.Properties.VariableNames) && ~(useNode && ismember('nearest_mesh_ind',tpoints.Properties.VariableNames))
                    
                    % Table of xyz
                    points = cat(2,tpoints.x, tpoints.y, tpoints.z);
                else
                    
                    % Table of mesh indices
                    assert(ismember('nearest_mesh_ind',tpoints.Properties.VariableNames), 'Table should have x column or nearest_mesh_ind');
                    mesh_ind = tpoints.nearest_mesh_ind;
                end
            elseif isrow(points) && numel(points) ~= 3
                % xyz
                points = points';
                
            elseif numel(points) > 1 && iscell(points) && any(size(points{1})==1)
                % cell array of mesh indices
                mesh_ind = points;
                
            end
            
            assert(~isempty(mesh_ind) || size(points,2) == 1 || size(points,2) == 3, 'plotPoint takes points with size n-by-1 or n-by-3'); 

            %^ ------ end argument handling ------ ^%
            
            %% get xyz
            if isempty(mesh_ind)
                if size(points, 2) == 1
                    surf = bp.get_psurf(surfName);
                    xyz = cat(1, surf.Vertices(points,:));
                elseif size(points, 2) == 3
                    xyz = points;
                end
            else
                % Map vertex sets to node
                lines = 0; % lines should be off when plotting vertex sets
                points = [];
                ndup = 0;
                for i = 1:numel(mesh_ind)
                    V = vector(mesh_ind{i});
                    surf = bp.get_psurf(surfName);
                    assert(~isempty(surf), 'Cannot find surface %s on current axis', surfName);
                    xyz = cat(1, surf.Vertices(V,:));
                    if mesh_clust_d > 0
                        clust = bp.clusterize(xyz, mesh_clust_d);
                        % for each cluster, take the xyz-mean and add it to points
                        for j = 1:max(clust)
                            points = cat(1, points, mean(xyz(clust == j, :), 1));
                            
                            % we are duplicating a point, so we should also duplicate other variables here....
                            if j > 1
                                % TODO
                                % duplicate any argument that is length n
                                fn = fieldnames(ip.Results);
                                for ifn = 1:numel(fn)
                                    param = fn{ifn};
                                    if ~ischar(ip.Results.(param)) && ~istable(ip.Results.(param)) && any(size(ip.Results.(param)) == n) && n > 1
                                        assert(exist(param, 'var')==1, 'trying to adjust parameter %s but it was not renamed exactly the same in function workspace', param);
                                        val = eval(param);
                                        val = vector(val);
                                        val = cat(1, val(1:i+ndup), val(i+ndup), val(i+1+ndup:end));
                                        eval([param '=val;']);
                                        
                                    end
                                end
                                ndup = ndup + 1;
                            end
                        end
                    else
                        % just usefirst index
                        points = cat(1, points, xyz(1,:));
                    end
                end
                xyz = points;
            end
            
            %% plot
            isRGB = isvector(color) && all(size(color) == [1 3]);
            clim = get(gca, 'CLim');
            if usePlot3
                [hsphere] = plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'o','MarkerFaceColor',color,'color',color,...
                    'MarkerSize',10,'linestyle','none', 'Tag','Electrode');
                hgroup = hggroup;
                hgroup.Parent = hsphere.Parent;
                hsphere.Parent = hgroup;
            else
                [hgroup, hsphere] = viz_view_sphere(xyz, 'color', color, 'radius', radius, 'opacity',alpha);
            end
            if isRGB
                % fix bug where matlab adjusts your color limit when plotting surfs, even if your surface is specified as RGB
                set(gca, 'CLim',clim);
            end
            if ~isempty(data) && numel(data) == 1
                for k = 1:length(hsphere)
                    hsphere.CData = repmat(data, size(hsphere.XData));
                    hsphere.FaceColor = 'flat';
                end
            end
            
            %% Line drawing
            if lines
                
                TOL = 0.2;
                LINE_WIDTH = 5;
                COL = color;
                d = pdist(xyz);
                
                % set default but check for conditions for neighbor method
                method = 'EXPECTED_DISTANCE'; 
                if ~isempty(element_info) && (istable(element_info) || exist(element_info, 'file'))
                    if isempty(tagName) && ~isempty(tpoints)
                        if ismember('tagName', tpoints.Properties.VariableNames) && numel(unique(tpoints.tagName)) == 1
                            tagName = unique(tpoints.tagName);
                        elseif ismember('chanName', tpoints.Properties.VariableNames) && numel(unique(util_split_stringnum(tpoints.chanName))) == 1
                            tagName = unique(util_split_stringnum(tpoints.chanName));
                        end
                    end
                    if ~isempty(tagName)
                        method = 'NEIGHBORS';                    
                    end
                end
                
                switch method
                    % Find all distances within TOL% distance of given lines distance (mm) and draw
                    % a line between those points
                    case 'EXPECTED_DISTANCE'
                        m = (d >= lines*(1-TOL) & d <= lines*(1+TOL));
                        m = squareform(m);
                        for i = 1:length(m)
                            for j=i+1:length(m)
                                if m(i,j) == 1
                                    l = line(xyz([i,j],1), xyz([i,j],2), xyz([i,j],3), 'color',COL,'lineWidth',LINE_WIDTH); 
                                    l.Parent = hgroup;
                                end
                            end
                        end
                    case 'NEIGHBORS'
                        if ~istable(element_info)
                            element_info = readtable(element_info);
                        end
                        neighbor_pairs = findNeighbors([],[], tagName, 1, 'element_info', element_info);
                        infoRow = element_info(strcmpi(element_info.tagName, tagName), :);
                        n_elecs = infoRow.bigDim*infoRow.smallDim;
                        elecs = setdiff(1:n_elecs, eval(infoRow.isCut{1}));
                        chanNames = util_combine_strnum(tagName, elecs);
                        
                        for r = 1:length(neighbor_pairs)
                            % these pairs are element number (not index). Find index:
                            % (assume xyz are in order)
                            i = neighbor_pairs(r,1);
                            j = neighbor_pairs(r,2);
                            i_name = sprintf('%s%d', tagName, i);
                            j_name = sprintf('%s%d', tagName, j);
                            i = find(strcmpi(chanNames, i_name));
                            j = find(strcmpi(chanNames, j_name));
                            
                            l = line(xyz([i,j],1), xyz([i,j],2), xyz([i,j],3), 'color',COL,'lineWidth',LINE_WIDTH); 
                            if isvalid(hgroup)
                                l.Parent = hgroup;
                            end
                        end
                        
                end
            end
            
            %% text-label, camera, legend
            if ischar(label) && ~iscellstr(label), label = {label}; end
            for i = 1 : length(label)
                cam = get(gca, 'CameraPosition');
                v = normalise(cam - xyz(i,:)); % vector from point to cam
                label_pos = xyz(i,:) + 2*v; % this pops label out a bit so it's visible
%                 t = text('Position',label_pos, 'String',label{i}, 'FontSize',13, 'HorizontalAlignment','center', 'color', labelColor);
%                 t.FontWeight = 'bold'; 
%                 t.Parent = hgroup;
%                 t.UserData = struct('root_xyz',xyz(i,:));
            end
            
            if ~isempty(legendLabel)
                hgroup.DisplayName = char(legendLabel);
            end
        end
        function plotPointGroup(bp, points, varargin)
            % plotPointGroup(bp, points, varargin) is similar to plotPoint, but it groups channels
            % by their hardware using colors and lines
            %
            % Inputs:
            %   points  - Either a list of xyz coords or table with x,y,z columns, or something else
            %               If you don't pass a table, make sure to pass chanNames or tagNames.
            %               See plotPoint for the full list of options!
            %
            % Optional (key-value) Inputs:
            %   element_info - if you want to draw the lines with the actual neighbors, pass subject's element info file path or table here
            %   chanNames   - cell of chanNames, 1 per point (may omit if points contains chanNames)
            %   tagNames    - cell of tagNames, 1 per point (may omit if points contains chanNames)
            %   ....        - Takes all optional arguments to plotPoint

            ip = inputParser;
            ip.addParameter('chanNames', []);
            ip.addParameter('tagNames', []);
            ip.KeepUnmatched = true;
            ip.parse(varargin{:});
            chans   = ip.Results.chanNames;
            tags    = ip.Results.tagNames;
            num     = [];
            
            % Table passed
            if istable(points)
                cols = points.Properties.VariableNames;
%                 if (ismember('x', cols) && ismember('y', cols) && ismember('z', cols))
%                     xyz = cat(2,points.x, points.y, points.z);
%                 elseif
%                     xyz = points;
%                 else
%                     error('Table passed without recognized columns');
%                 end
                
                if ismember('chanName', cols)
                    chans = points.chanName;
                    [tags,num] = util_split_stringnum(chans);
                    
                elseif ismember('tagName', cols)
                    chans = [];
                    tags = points.tagName;
                end
            end
            
            utags = unique(tags, 'stable');
            n = numel(utags);
            rng(bp.rngseed);
%             colors = hsv(n);
%             colors = colors(randperm(n),:);

            colors = parula(n); 

            for i = 1:n
                ndx = strcmpi(tags, utags{i});
                newargs = {'color',colors(i,:), 'legend',utags{i}, 'lines',10, 'tagName', utags{i}};
                
                if ~isempty(num)
                    newargs = [newargs {'label',arrayfun(@(x) {x}, num(ndx))}];
                end
                
                newargs = [newargs varargin];
                bp.plotPoint(points(ndx,:), newargs{:});
            end
            
            if ~ismember('element_info', fieldnames(ip.Results))
                bp.hint('If the connected lines look wrong, try passing element_info so the plotter knows how to actually connect them!')
            end
        end
        function region_surf = plotRegion(bp, V, varargin)
            % plotRegion(bp, V, varargin)
            %
            % Colors the set of vertices given by indices in V
            %
            %  Input:
            %   V       - This can be a list of node indices that define the region you are plotting.
            %
            % Optional key-Value Pairs:
            %   surf    - Name of surface on which to plot
            %   color   - rgb color of region OR data value (indexed color). single value or 1 row per vertex
            %   blend   - opacity in [0,1]. By default this is 0.25 if ~useNewSurf and 1.00 if useNewSurf
            %   boundaryWidth - pass a boundary width. Default no boundary. (2 is good)
            %   useNewSurf - true to create a new surface for the region
            %   haxis   - axis handle
            %
            % The bp surface will get colored at all points in V. useNewSurf will
            % instead leave bp.surfaces unchanged and overlay a translucent new surface.
            % UseBoundary always creates a new surface with just the boundary colored.
            % If you use boundary, appearance is better if you also useNewSurf,
            % although you will be adding 2 new patches which will slow performance if
            % you are plotting a lot.
            %
            %
            % More details on surface colors:
            %   Colors are determined by the surface (patch object) property FaceVertexCData. This must
            %   be either a single RGB value or (more likely) an N x 3 or N x 1 matrix of values, where
            %   N is the number of vertices (198,812 for std mesh). The number of columns determines color type, which
            %   can be either of the following two modes:
            %
            %       True Color - each vertex gets an RGB triplet. This is the only way to truly blend colors if you
            %                    want to implement your own blending for overlays (e.g. mean(base color, overlay color)).
            %                   The disadvantage is that RGB is static in relation to the colormap
            %       Indexed color - each vertex gets a single data value, which indexes into the axis's current colormap.
            %                   You cannot mix two colors or two colormaps, but the color changes with the color map.
            % 
            % Layering surfaces:
            %       To get around the layering limitation of indexed colors, you could create layers by creating a new surface and adjust
            %       the transparency of it. In this way, the matlab renderer will calculate the blend for you. The new surface is still 
            %       limited by the above problems though. For example, you still can't have two different color maps on the same axis.

            
            region_surf = [];

            if ~isnumeric(V)
                keyboard;
            end
            
            %% --- Argument handling --- %%
            ip = inputParser;
            ip.addParameter('surf',[]);
            ip.addParameter('color',[0 0 1]);
            ip.addParameter('blend',[]);
            ip.addParameter('name','boundary');
            ip.addParameter('boundaryWidth',0);
            ip.addParameter('useNewSurf',0);
            ip.addParameter('haxis',[]);
            ip.parse(varargin{:});
            surfName    = ip.Results.surf;
            color       = ip.Results.color;
            blend       = ip.Results.blend;
            name        = ip.Results.name;
            boundaryWidth = ip.Results.boundaryWidth;
            useNewSurf  = ip.Results.useNewSurf;
            haxis       = ip.Results.haxis;
    
            if isempty(haxis)
                haxis = gca();
            end
            
            
            surf = bp.get_psurf(surfName, haxis);
            assert(~isempty(surf), 'plotRegion::surf is empty')

            if isempty(blend)
                % default blends
                if useNewSurf,  blend = 1;
                else,           blend = 0.25;
                end
            end         
            %% --- Argument handling end --- ^%
            
            % Add color
            V = unique(V, 'stable');
                        
            if size(color, 2) == 3
                % rgb - blend with existing rgb
                fvcd = surf.FaceVertexCData;
                if size(fvcd, 2) == 1
                    % convert current indexed-color to rgb
                    bp.saveToRGB(surf);
                    fvcd = surf.FaceVertexCData;
                end
                
                if size(color,1) == 1
                    fvcd(V, :) = blend * repmat(double(color), length(V), 1) + fvcd(V,:) * (1-blend);
                else
                    fvcd(V, :) = blend * double(color) + fvcd(V,:) * (1-blend);
                end
                % if ~useNewSurf
                %%%%%%
                % JD editing here.
                %%%%%%
                % Not really understanding this code, but from what I do
                % understand, unless we express this equality, our nx3
                % color data will never be used. 
                    surf.FaceVertexCData = fvcd;
                % end
                use_rgb = 1;
            
            elseif size(color,2)==1 && size(color,1)==length(V) && ~useNewSurf
                % non-RGB
                surf.FaceVertexCData(V) = color;
                
                use_rgb = 0;
                
            else
                use_rgb = 0;
            end

            
            if useNewSurf || boundaryWidth > 0
                
                
                
                if boundaryWidth > 0
                    % Only color boundary
                    [Vb,Fb] = bp.get_region_boundary(V, surf);
                    rgb = nan(length(surf.Vertices), 3);
                    rgb(Vb,:) = repmat(color,length(Vb),1);
                    
                    bnd_surf = patch(haxis, ...,
                        'faces',int32(Fb),...
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
                        CData = surf.FaceVertexCData;
                        % CData(V,:) = repmat(color, length(V),1);
                    else
                        if isscalar(color) || bp.stdNumVert == length(color)
                            CData = color;
                        else
                            CData = zeros(bp.stdNumVert, 1);
                            CData(V) = color;
                        end
                    end
                    
                    % Subsample faces and reindex its indices
                     fmask = ismember(surf.Faces,V);
                     f = surf.Faces(sum(fmask,2) == 3, :);
%                     F = indexOf(f,V);


                    region_surf = patch(haxis, ...,
                        'faces',f,...
                        'vertices',surf.Vertices,...
                        'facevertexcdata',full(CData),...
                        'edgecolor','none',...
                        'facecolor','interp',...
                        'faceAlpha', blend, ...
                        'specularstrength',0, ...
                        'faceLighting', 'gouraud', ...
                        'tag','region');
                    
                else
                    region_surf = [];
                end
                
            end
        end
        function region_surf = plotdset(bp, data, varargin)
            % This wrapper around plotRegion plots the data values given in Data
            % Data must be the same length as the number of vertices
                region_surf = bp.plotRegion(1:length(data), 'color',data, varargin{:});

            
        end
        function region_surf = plotAtlas(bp, fname, varargin)
            % plotAtlas(bp, fname, varargin) colors a surface according to an atlas file
            %
            %   This is a wrapper around plotdset (which calls plotRegion so all arguments get passed to plotRegion)
            %
            % Example: bp.plotAtlas('SUMA/std.141.lh.aparc.a2009s.annot.niml.dset', 'surf','lh')
            %
            % Hint - Use the data cursor to display the labels wherever you point your mouse!
            %
            % Note: we only have SUMA data for destrieux ("a2009s") atlas currently. The other atlases are not
            % yet ported to SUMA from Freesurfer
            
            
            temp = which('afni_niml_annotation');
            assert(~isempty(temp), 'Need the afni library for reading afni_niml_annotation');
            try
                [~,label,ctab] = afni_niml_annotation(fname);
                [~,rgb_ndx] = ismember(label, ctab.table(:,5));
            catch
                x = afni_niml_read(fname);
                label = x{1}.nodes{1}.data;
                keyboard;
                
            end
            rgb_ndx(rgb_ndx == 0) = 1;
            names = ctab.struct_names(rgb_ndx);
            rgb = ctab.table(rgb_ndx, 1:3);
            rgb = double(rgb) / 255;
            region_surf = bp.plotdset(rgb, varargin{:});
            
            % Set up datacurser
            f = gcf();
            dcm = datacursormode(f);
            dcm.UpdateFcn = [];
            dcm.UpdateFcn = {@bp.dataCurserAtlas, names};
        end
        
        function plotSulc(bp, bd, surfName)
            % plotCurvature(bp, bd, surfName) Plots a braindata object's curvature on the surface
            %
            %
            % Input
            %   bd      - braindata2 object OR path the a directory containing surfaces. Default pwd
            %   surf    - Name of loaded surface. Optional if only 1 surface is loaded
            %
            % Details:
            %   the curvature is "snapped" to -1/+1 using heaviside
            
            if nargin < 2 || isempty(bd), bd = pwd; end
            if nargin < 3, surfName = ''; end
            
            if isa(bd, 'braindata2')
                dir_surfaces = bd.filepaths.surfaces;
            else
                dir_surfaces = bd;
            end
            
            surf = bp.get_psurf(surfName);
            
            % Set some defaults
            colormap('gray')
            set(gca,'clim',[-1,1.5]); %(this ensures that white will actually be grayish)
            
            if strfound(surf.Tag, 'lh')
                hemi = 'lh';
            elseif strfound(surf.Tag, 'rh')
                hemi = 'rh';
            else
                error('No hemisphere is surface tag -- I dont know what file to load');
            end
            
            
            fname = fullfile(dir_surfaces, sprintf('std.141.%s.sulc.niml.dset', hemi));
            data = bp.dsetToRegion(fname);
            data = 2*heaviside(-data) - 1;
            % note that we negate the sulc values because we want sulcus to be dark
            
            bp.plotdset(data, 'surf',surfName);
            
            % Data->RGB (because this is a base layer and shouldn't be colormap-dependent)
            bp.saveToRGB(surf);
            
        end
        
        function [vertexData, vertexCnt, VPerROIMat, region_surf] = plotRegionsDataMat(bp, VPerROIMat, valPerROIC, VFind, varargin)
            % plotRegionsDataMat plots ROIC data quickly
            % 
            % USAGE
            %
            % INPUTS
            %   VPerROIMat  - This can be either a VPerROIMat OR a VPerROI cell array
            %   valPerROIC - This is a sparse matrix array of length nROIC (or 2xnROIC). See braindata2/leadData2ROICMat.
            %   V
            %
            % OPTIONAL INPUT
            %   surf        - Surface name to plot on
            % 
            % OUTPUTS
            %   vertexData  - sparse array of data plotted at each vertex
            %   vertexCnt   - how many data points are at each vertex
            %   VPerROIMat  - VPerROIMat, sparse matrix
            %   region_surf - handle to patch of plotted region
            %
            % DETAILS
            %
            
            % input
            ip = inputParser;
            ip.addParameter('surf',[]);
            ip.addParameter('haxis',[]);
            ip.parse(varargin{:});
            
            surfName = ip.Results.surf;
            
            % output
            vertexCnt = [];
            
            assert(isvector(valPerROIC), 'valPerROICMat must be a vector (ie 1 value per ROIC), not a matrix');
            valPerROIC = valPerROIC(:);
            assert( size(VPerROIMat,2) == numel(valPerROIC), 'VPerROIMat columns (nROI, %d) must match valPerROICMat length (%d)', ...
                    size(VPerROIMat,2), numel(valPerROIC));
            
            
            

            % Multiply (ie ROI-data average at each vertex)
            vertexData = VPerROIMat * valPerROIC;
            %vertexData = VPerROIMat(:,sparseMask) * valPerROIC(sparseMask);
                        
            % Plot
            region_surf = bp.plotRegion(VFind, 'color',vertexData, 'surf',surfName, 'useNewSurf',1, varargin{:});
            
        end
        
        function [vertexData, vertexCnt] = plotRegionsData(bp, VPerROI, valPerROIC, varargin)
            % [vertexData, vertexCnt] = plotRegionsData(bp, VPerROI, valPerROIC, varargin)
            % plots data on multiple vertex sets
            %
            % Input:
            %   VPerROI     - n x 1 cell array of surface index arrays; Vertices for each ROI
            %   valPerROIC  - n x 1 array of data values to plot; Value for each ROI center
            %
            % Optional Key-value pairs:
            %   surf        - Surface name to plot on
            %   rh_begin    - If you're giving [lh; rh], begin rh at this index
            %   clim        - [low high]
            %   cmap        - colormap
            %   blend       - [default=1] 1 gives full color, 0.1 gives 10 percent color
            %   cumulative  - if 0 [default], take mean of overlapping vertex sets
            %                 if 1, take sum of overlapping vertex sets
            %                 if 2, take max
            %                 if 3, median
            % Output:
            %   vertexData  - data plotted at each vertex from [1:max(vertexSets)]
            %   vertexCnt   - how many data points are at each vertex
            %
            % Details:
            %   Uses plotRegion to plot a new surface and uses colormap-indexed colors
            %
            %
            
            
            ip = inputParser;
            ip.addParameter('surf',[]);
            ip.addParameter('rh_begin',[]);
            ip.addParameter('clim',get(gca,'CLim'));
            ip.addParameter('cmap',colormap);
            ip.addParameter('blend',1);
            ip.addParameter('cumulative',0);
            ip.addParameter('haxis', []);
            ip.parse(varargin{:});
            
            surfName = ip.Results.surf;
            rh_begin = ip.Results.rh_begin;
            clim = ip.Results.clim;
            cmap = ip.Results.cmap;
            blend = ip.Results.blend;
            cumulative = ip.Results.cumulative;
            haxis = ip.Results.haxis;
            if isempty(haxis)
                haxis = gca();
            end
            % Check for lh/rh split
            
            
            %%%%%
            % JD edits: 
            %%%%%%
            % if rh_begin > 1
            if rh_begin > 0 
                modargin = [varargin {'surf','lh','rh_begin',0}];
                [vertexData_l, vertexCnt_l] = bp.plotRegionsData(VPerROI(1:rh_begin-1), valPerROIC(1:rh_begin-1), modargin{:});
                modargin = [varargin {'surf','rh','rh_begin',0}];
                [vertexData_r, vertexCnt_r] = bp.plotRegionsData(VPerROI(rh_begin:end), valPerROIC(rh_begin:end), modargin{:});
                vertexData = cat(1, vertexData_l, vertexData_r);
                vertexCnt = cat(1, vertexCnt_l, vertexCnt_r);
                return;
            elseif isempty(surfName)
                %%%% 
                % JD edits
                %%%%
                return
            end
            
            % main
            n = length(VPerROI);
            assert(length(valPerROIC) == n, 'data and VPerROI must be the same size');
            %surf = bp.get_psurf(surfName);
            

            vertexData = zeros(bp.stdNumVert, 1);
            vertexCnt  = zeros(bp.stdNumVert, 1);
            vertexDataC= cell(bp.stdNumVert, 1);
            
            
            % aggregate data from different vertex sets onto each vertex
            for i = 1:n
                vertexCnt(VPerROI{i}) = vertexCnt(VPerROI{i}) + 1;
                if cumulative == 2
                    % Max
                    vertexData(VPerROI{i}) = nanmax(vertexData(VPerROI{i}), valPerROIC(i));
                elseif cumulative == 3
                    % Cells
                    error('Not implemented: median (cumulative==3)')
                    %vertexDataC(VPerROI{i}) = [vertexDataC(VPerROI{i}) repmat({vertexData(VPerROI{i})}, numel(VPerROI{i},1)];
                else    
                    % Sum (for mean, divide later)
                    vertexData(VPerROI{i}) = vertexData(VPerROI{i}) + valPerROIC(i);
                end
            end
            
            % average at each vertex
            findMask = vertexCnt > 0;
            vertexData(~findMask) = nan;
            if cumulative == 0 % average the sum
                vertexData(findMask) = vertexData(findMask) ./ vertexCnt(findMask);
            elseif cumulative == 3
                %vertexData(findMask) = cellfun(@nanmedian, 
                error('Not implemented')
            end
            
            
            %%%%%%%%
            % JD edits here. Adding support for clim
            if ismember('clim',varargin(1:2:end))
                
                vertexData = vertexData(findMask);
                % clim has been passed.
                % We'll have to convert our colors into RGB. 
                
                % Put it all in bounds 
                vertexData = max(clim(1),min(vertexData,clim(2)));
                % Normalize to bounds
                vertexData = (vertexData - clim(1)) / range(clim);
                % Times 64
                vertexData = round(vertexData * size(cmap,1));
                
                vertexData = squeeze(ind2rgb(vertexData,cmap));
            end
            %%%%%%%%%
                
            
            
            bp.plotRegion(find(findMask),'useNewSurf',1,'color',vertexData, 'blend',blend, 'surf',surfName, 'haxis',haxis);
        end
        function replotLabels(bp)
            % Takes all text labels and makes them visible by repositioning them between the camera and their 
            % root location
            txts = findall(gca, 'type', 'text');
            
            POP_OUT_MM = 5;
            
            for i = 1:numel(txts)
                txt = txts(i);

                % Note, if the texts were plotted with plotPoint, they should have the root xyz in userData 
                if ~isempty(txt.UserData) && isfield(txt.UserData, 'root_xyz')
                    xyz = txt.UserData.root_xyz;
                    cam = get(gca, 'CameraPosition');
                    v = normalise(cam - xyz); % vector from point to cam
                    label_pos = xyz + POP_OUT_MM * v; % this pops label out a bit so it's visible
                    txt.Position = label_pos;
                end
            end
            
            
            
        end
    end
    
    methods % Viewing
        
        
        function varargout = viewAll(bp, fig, copyFromAx)
            % viewAll(bp, [fig, copyFromAx]) takes the current axis and duplicates it onto another figure in standard views
            % Useful for creating a figure to save
            if nargin < 3, copyFromAx = gca(); end
            if nargin < 2 || isempty(fig), fig = figure(); end

            RMV_LIGHTS = 1; % remove the ensemble of lights, but put one direct light at camera
            
            clf(fig);
            figure(fig);
            pause(1);
            fig.Units = 'Normalized';
            views = {'ventral', 'anterior', 'dorsal', 'left', 'right'};
            for i=1:5
                nax = subplot(2,3,i);
                copyobj(copyFromAx.Children, nax);
                bp.apply_axis_settings();
                bp.view(views{i});
                bp.replotLabels();
                
                if RMV_LIGHTS
                    lights = findall(nax, 'type', 'light');
                    for j=1:numel(lights)
                        lights(j).Parent = [];
                    end
                    l = camlight('infinite');
                    
                end
                axList(i) = nax;
            end

            % Axis location adjustments
            axList(1).Position([2,4]) = [0.5, 0.5];
            axList(3).Position([2,4]) = [0.5, 0.5];
            axList(2).Position([2,4]) = [0.3, 0.5];
            axList(4).Position([1,3]) = [0, 0.50];  
            axList(5).Position([1,3]) = [0.5, 0.50];
            
            if ~isempty(legend())
                bp.legend('orientation', 'horizontal');
            end
            
            if nargout >= 1
                varargout{1} = fig;
            end
        end
        
    end
    
    methods (Static) % Viewing
        function view(name)
            % view(bp, name) sets the view to a given name
            % NAME can be one of:
            %   left   (anterior left)
            %   right  (anterior right)
            %   dorsal          (anterior top)
            %   ventral         (anterior top)
            %   anterior        (anterior top, LH = right)
            %   
            if nargin < 1, name = 'dorsal'; end
            names = {'left' 'right' 'dorsal', 'ventral', 'anterior'};
            
            assert(ismember(lower(name), names), 'Incorrect vew name, see help');
            switch lower(name)
                case 'left',        view(-90,0);
                case 'right',       view(90,0);
                case 'dorsal',      view(0,90);
                case 'ventral',     view(180,-90);
                case 'anterior',    view(180,0);
            end
            
        end
        function showAxis()
            % Displays axis box with labels
            axis on;
            xlabel x;
            ylabel y;
            zlabel z;
        end
        function rotate(duration, T, fps, useVideo, varargin)
            if nargin < 1 || isempty(duration)
                duration = 10; % seconds
            end
            if nargin < 2 || isempty(T)
                T = 4; % seconds per 2pi (full) revolution
            end
            if nargin < 3 || isempty(fps)
                fps = 60;
            end
            if nargin < 4 || isempty(useVideo)
                useVideo = 0;
            end
            warning('off');
            
            assert(~useVideo, 'Not supported yet');
            
            axis vis3d;
            ax = gca();
            tstep = 1/fps;
            
            t0 = tic();
            t = toc(t0);
            while t < duration
                view(360 * t/T, 0);
                pause(tstep);
                t = toc(t0);
                
                if contains(varargin,'disco')
                    brainplotter.camlights(8,'easter');
                end
            end
            warning('on');
        end
        function rgb = saveToRGB(p)
            % This function RGB-ifies the given patch's colordata
            % See brainplotter/plotRegion for more details on patch colors
            % This is used when you plotted a patch as an indexed color but want to convert it to RGB
            %
            % INPUT
            %   p - matlab patch that you want to RGB-ify
            %
            % See Also: plotRegion

            cdata = p.FaceVertexCData;
            if size(cdata, 2) ~= 1
                return;
            end
            mask_nans = isnan(cdata);
            
            % Get the current values in terms of the current colormap
            % the getcc functions map data into 1:n
            cm = colormap();
            clim = get(gca, 'CLim');
            n = size(cm,1);
            getc = @(val) round( (n-1) * (val-clim(1)) / (clim(2)-clim(1)) + 1 );
            getcc = @(val) min(max(getc(val), 1), n);
            ndx = arrayfun(getcc, cdata);
            rgb = cm(ndx, 1:3);
            rgb(mask_nans,:) = nan(sum(mask_nans), 3);
            
            p.FaceVertexCData = rgb;
        end
    end
    
    methods % Theoretically private
        
        
        
        
    end
    
    methods (Static)
        function version()
            fprintf('Version: Brainplotter 08/2018\n');
        end
        function [vertexData, sparseMask] = dsetToRegion(path2dset)
            % dsetToRegion(path2dset) maps an afni dset onto a surface
            data = afni_niml_readsimple(path2dset);
            vertexData = data.data;
            sparseMask = ones(size(vertexData));
        end
        function lights = camlights(intensity, varargin)
            % lights = camlights(intensity)
            if nargin < 1 || isempty(intensity), intensity = 5; end
            
            
            lights = findall(gca,'type','light');
            for i=1:length(lights), lights(i).Parent=[]; end
            
            lights = {};
            for el = [45, -35]
                for az = [0, 90, 190, 270]
                    lights = [lights {lightangle(az,el)}];
                end
            end
            
            cols = [0 1 1; 0 0 1; 1 0 1; 0 1 0; [1 0 0]];
            
            for i=1:length(lights)
                
                if any(strcmpi(varargin, 'easter'))
                    lights{i}.Color = cols(randi(length(cols)), 1:3) .* intensity ./ length(lights);
                else
                    lights{i}.Color = [1 1 1] ./ length(lights) .* intensity;
                    lights{i}.Style = 'infinite';
                end
            end
            if contains(varargin, 'disco')
                for i=1:1000, brainplotter.camlights(5,'easter'); pause(0.1); end
            end
        end
        function [psurf, name] = get_psurf(name, haxis)
            % psurf = get_psurf(name)
            % selects all plotted patches with any part of their tag matching the given name
            
            if nargin < 1, name = []; end
            if nargin < 2, haxis = gca(); end
            
            % remove duplicates surfaces (not boundaries)
            psurf = findall(haxis,'type','patch');
            if ~isempty(psurf)
                tags = psurf.get('Tag');
                if ischar(tags)
                    tags=cellstr(tags); 
                end
                tags = tags(~strcmpi(tags,'boundary'));
                [~,i_unique] = unique(tags,'stable');
                i_dup = setdiff(1:numel(tags), i_unique);
                for i = i_dup
                    if ~isempty(tags{i})
                        psurf(i).Parent = [];
                    end
                end
            end
            
            if isempty(name)
                psurf = findall(haxis,'type','patch'); 
            else
                
                psurf = findall(haxis, '-regexp', 'Tag', name);
                if numel(psurf) > 1
                	psurf = psurf(strfound(psurf.get('Tag'), name));
                elseif numel(psurf) == 0
                    %keyboard;
                end
            end
            
            
                
            if numel(psurf) > 1
                brainplotter.hint('There are >1 plotted patches on the axis; I do not know which to plo1 (try passing a param like: ''surf'',''pial_rh'')');
                i = menuText('Which one?', {psurf.Tag}, 'multiSelect',0);
                psurf = psurf(i);
            elseif numel(psurf) == 0
                return;
            end

            name = psurf.Tag;
        end
        function setOpacity(opacity)
            % setOpacity(opacity) sets alpha value of surface to opacity
            % opacity is in [0,1] with 1 being opaque

            ax = gca;
            patches = findobj(ax, 'Type', 'patch');
            for p = 1 : length(patches)
                patches(p).FaceAlpha = opacity;
            end
            
        end
        function setEdgeColor(edgecolor)
            % setOpacity(opacity) sets alpha value of surface to opacity
            % opacity is in [0,1] with 1 being opaque
            
            ax = gca;
            patches = findobj(ax, 'Type', 'patch');
            for p = 1 : length(patches)
                patches(p).EdgeColor = edgecolor;
            end
            
        end
        function setFaceColor(facecolor)
            % setOpacity(opacity) sets alpha value of surface to opacity
            % opacity is in [0,1] with 1 being opaque
            
            ax = gca;
            patches = findobj(ax, 'Type', 'patch');
            for p = 1 : length(patches)
                patches(p).FaceColor = facecolor;
            end
            
        end
        
        
        function clearPoints()
            % clearPoints removes all plotted electrodes from the current axis by removing:
            %   Text
            %   Lines
            %   Objects tagged with "electrode"
            points = findall(gca, '-regexp', 'Tag','electrode');
            for i=1:length(points)
                points(i).Parent = [];
            end
            lines = findall(gca, 'type', 'Line');
            for i =1:length(lines)
                lines(i).Parent = [];
            end
            
            text = findall(gca, 'type', 'text');
            for i =1:length(text)
                text(i).Parent = [];
            end
            
            
        end
        function varargout = legend(varargin)
            %   legend([...]) puts a legend on the current axis. It is a wrapper around the normal legend function
            %
            %   The changes that this wrapper makes are as follows:
            %       - Only axis objects with a non-empty DisplayName property are included
            %       - A height and width are set (according to the given orientation)
            %       - Latex interpreter is off
            %       - No duplicates of objects with the same DisplayName
            %
            %   See also  LEGEND


            DEFAULT_HORIZONTAL_HT = 0.04; % use lgd.Height to change
            default_args = {'orientation','vertical', 'Location','best'};
            newargin = [default_args varargin];
            
            hAxes = gca;
            
            % Look for all "group" objects (calls to plotPoint are put into groups)
            % Filter those groups to only those groups with a DisplayName property set
            hgroups = findobj(hAxes, 'type', 'hggroup'); 
            hgroups = hgroups(~cellfun('isempty', {hgroups.DisplayName}));
            if isempty(hgroups)
                return
            end
            
            seen_names = [];
            clear legendObj legendObjs
            for igroup = 1:length(hgroups)
                displayName = get(hgroups(igroup), 'DisplayName');
                if ~ismember(displayName, seen_names)
                    legendObj            = findall(hgroups(igroup), 'type', 'surface');
                    if isempty(legendObj)
                        legendObj        = findall(hgroups(igroup), 'type', 'line');
                    end
                    legendObj            = legendObj(1);
                    
                    legendObjs(numel(seen_names)+1) = legendObj;
                    seen_names{numel(seen_names)+1} = displayName;
                end
                
            end
            hLegend = legend(legendObjs, seen_names{:},  newargin{:});
            
            hLegend.Interpreter = 'none';
            if strcmpi(hLegend.Orientation, 'horizontal')    
                hLegend.Position = [0, 0, 1, DEFAULT_HORIZONTAL_HT];
            end
            
            if nargout == 1
                % Return array of handles to all legends created if asked for
                varargout = hLegend;
            end
        end
        
        function d = dist_geodesic(gsurf, ndx_seed, dest_ndx)
            % Computes geodesic distance from one surface point to all other points on the surface
            %
            % Inputs:
            %   gsurf - Gifti surface object (with vertices and faces fields)
            %   ndx_seed - Index into surface mesh
            %
            % Optional Input:
            %   dest_ndx - if you want the distance to a specific node, pass that node here (faster)
            %
            % Example:
            %   d = brainplotter.dist_geodesic(bp.surfaces.pial_lh, mesh_ndx)
            %   d = brainplotter.dist_geodesic(bp.surfaces.pial_lh, mesh_ndx1, mesh_ndx2)

            if nargin < 3
                dest_ndx = [];
            end
            
            V = double(gsurf.vertices)';
            F = double(gsurf.faces)';
            
            if isempty(dest_ndx);
                d = surfing_dijkstradist(V, F, ndx_seed);
            else
                n2f = surfing_nodeidxs2faceidxs_custom(gsurf.faces', gsurf.vertices');

                d_euclid = norm(double(gsurf.vertices([ndx_seed; dest_ndx],:)));
                PADCLIP = 1.10;
                % The following clips the surface using euclidean distance (with a 10 percent padding) for speed improvement
                % and then computes geodesic distance with fast-marching
                
                
                [mesh_ndxs, dists, ~, ~] = surfing_circleROI_custom(V, F, ndx_seed, d_euclid*PADCLIP, 'geodesic',n2f);
                
                %%%%%
                % MT version
                % d = dists(mesh_ndxs == dest_ndx);
                %%%%%%
                
                %%%%%%
                % JD Version
                d = dists(ismember(mesh_ndxs, dest_ndx));
                % Alternatively:d = dists(any(mesh_ndxs == dest_ndx));
                
                % For some reasons currently unknown to me, the above call
                % can fail, that is, elements of dest_ndx exist which
                % are not part of mesh_ndxs. Have yet to figure it out but
                % as of now we'll provide a patch hack. 
                                
                if length(d) < length(dest_ndx)
                    % d = ones(size(dest_ndx)) * -99;
                    d = nan(size(dest_ndx)); 
                end
                %%%%%%%%%%%
            end
            
        end
    end
    
    methods (Static) % And theoretically private
        function [V,F] = get_region_boundary(V, surf)
            % V is index into surf.vertices
            % F is index into surf.faces
            
            
            F = double(surf.Faces);
            
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
        function apply_surf_settings()
            % Lighting and surfaces
            % These calls affect all surfaces, so if called after electrodes are plotted, electrodes get recolored
            % material dull;
            material([0.3 0.6 0.5 10 0.5])
            shading interp;
            lighting phong; % phong, flat 
            brainplotter.camlights();
        end
        function apply_axis_settings()
            % standard axis / view settings
            
            axis tight equal off;
            set(gcf, 'color', 'w') 
            %set(gcf, 'color', 'k') 
            set(gca, 'Clipping', 'off');
            
            xlabel x;
            ylabel y;
            zlabel z;
            hold on;
        end
        function hint(txt)
            % Display a hint to user in a standard way
            disp('');
            brainplotter.line();
            ll = 95; % line length
            n = length(txt);
            nlines = ceil(n/ll);
            msg = '';
            for i=1:nlines
                if i==1,                msg = [msg, '\t'   txt(ll*(i-1)+1 : min(ll*i,n))]; 
                elseif  i < nlines,     msg = [msg, '-\n\t\t-' strtrim(txt(ll*(i-1)+1 : ll*i)), '\n']; 
                else,                   msg = [msg, '-\n\t\t-' strtrim(txt(ll*(i-1)+1 : ll*(i-1) + mod(n,ll))), '\n']; 
                end
            end
            fprintf(['  ** Hint! ** ', msg '\n']);
            brainplotter.line();
            disp('');
        end   
        function line(title)
            if nargin < 1, title=0; end
            if title
                fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Brain Plotter ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
            else
                fprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
            end
        end
        function clust = clusterize(xyz, d)
            % clusters xyz points into clusters with centers farther than (1+TOL)*d away
            % clust will be integer array of length matching xyz with integers in 1:k for k clusters
            %
            % e.g. for electrode with r=1.5 mm, use 3 for d
            % 
            % Algorithm (by Mike)
            %   all points within (1+TOL)*d from a point are in that point's cluster
            %
            TOL = 0.1;
            if nargin < 2, d = 3; end
            
            dmax = (1+TOL)*d;
            n = length(xyz);
            D = pdist(xyz);
            M = squareform(D);
            
            clust = zeros(n,1);
            front = zeros(n,1);

            % Note: single letters are indices, words are masks
            k = 0;
            
            while ~all(clust)
                % Add an unclustered element to the frontier
                k = k+1;
                i = find(~clust,1);
                clust(i) = k;
                front(i) = 1;
                
                while any(front)
                    
                    % get any element in the frontier and remove it from the frontier
                    f = find(front, 1);
                    front(f) = 0;
                   
                    % add all unassigned points closer than the maximum to this cluster
                    % and add those points to the frontier
                    toAdd = vector((M(f,:) <= dmax)) & ~clust;
                    clust(toAdd) = clust(f);
                    front(toAdd) = 1;
                end
            end
            
            
        end
        function txt = dataCurserAtlas(h, event, names)
            [~,ndx] = min(sqrt(sum((event.Target.Vertices - event.Position) .^ 2, 2)));
            %txt = atlasFSTranslate('destrieux', names{ndx(1)}, 'translation','pretty');
            txt = atlasFSTranslate('desikan', names{ndx(1)}, 'translation','short');
        end
    end
end
    
    
