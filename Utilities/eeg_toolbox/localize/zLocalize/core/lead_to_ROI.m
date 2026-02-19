function lead_roi_lut = lead_to_ROI(roi_lut, proj_coords, pial_surf, envelope_surf, roi_radius, varargin)
% lead_to_ROI 
%
% Usage: lead_to_ROI = lead_to_ROI(roi_lut, proj_coords, pial_surf, envelope_surf, roi_radius, working_dir, ...)
%
% Note: Make sure that you are passing in the hemisphere surface corresponding to the mesh indices
%       (remember, mesh indices are not symmetric between hemispheres)
%
% Inputs:
%   roi_lut         - output from grow_ROIs: ROI table with columns: ROIC_mesh_ndx, vertex, d
%   proj_coords     - output from proj_coords: Table with xyz coordiantes of every electrode
%   pial_surf       - A struct with faces and vertices fields representing the pial surface mesh
%   evenlope_surf   - A struct with faces and vertices fields representing the cortical envelope mesh
%   roi_radius      - Radius of ROIs (geodesic, mm)
%   
% Optional Input:
%   working_dir     - directory where intermediate files and output files will be saved
%   lead_mesh_lut   - table with chanName, nearest_mesh_ind, geodisc_mesh_ind, geodisc_dist columns
%
% Input key-values:
%   log_file        - path to log file. 
%   contact_radius  - radius to model electrode as disk on cortical envelope. Default 1.5 (mm).
%
% Output: 
%   lead_roi_lut    - table of electrodes' ROIs. This table will be of the same form as
%                     the mr_coords table with new columns: lead_verts, ROIC_mesh_ndx, which are arrays of
%                     pial surface mesh indices.
%
% File Output:
%   
% Description:
%   Takes XYZ coordinates of electrodes and the subject brain/ROIs and returns the ROIs to 
%   which each electrode maps
%   
%

% Revision History:
%   03/17 MST   - Created
%   01/18 MST   - envelope_surf not used. Add lead_mesh_lut

%  Copyright (C) 2017 Mike Trotta

    % -------------
    %%--- Setup ---
    % -------------
    
    % Input
    ip = inputParser;
    ip.addOptional('working_dir', pwd, @ischar);
    ip.addOptional('lead_mesh_lut', [], @istable);
    ip.addParameter('log_file', []);
    ip.addParameter('contact_radius', 1.5);
    ip.parse(varargin{:});
    log_file        = ip.Results.log_file;
    contact_radius  = ip.Results.contact_radius;
    lead_mesh_lut   = ip.Results.lead_mesh_lut;
    
    % ---------------
    % Check input
    % ---------------
    if ismatrix(proj_coords) && size(proj_coords,2)==3
        proj_coords = array2table(proj_coords, 'VariableNames', {'x','y','z'});
    end
    
    required_cols = {'x','y','z'};
    assert(istable(proj_coords) && length(required_cols) == length(intersect(required_cols, proj_coords.Properties.VariableNames)),...
        'proj_coords must be a table with x, y, and z columns');

    required_cols = {'vertex','ROIC_mesh_ndx','d'};
    assert(istable(roi_lut) && length(required_cols) == length(intersect(required_cols, roi_lut.Properties.VariableNames)),...
        'roi_lut must be a table with vertex, ROIC_mesh_ndx, and d columns');
    
    required_fields = {'faces','vertices'};
    if sum(ismember(required_fields, fieldnames(pial_surf))) < 2 ...
        error('Pial and evelope surfaces must be structures with .faces and .vertices fields');
    end
    
    d_max = max(roi_lut.d);
    if d_max + 1 < roi_radius
        warning('Given ROI radius (%d) is greater than greatest distance in roi_lut (%d) (radius should match that from grow_ROIs)',...
            roi_radius, d_max);
    end
    
    pial_surf.vertices = single(pial_surf.vertices);
%     if ~isempty(envelope_surf)
%         envelope_surf.vertices = single(pial_surf.envelope_surf);
%     end
    
    
    % ---------------
    %%--- Mapping ---
    % ---------------
    
    % Map electrode to pial nodes (via sphere model)
    lead_roi_lut = proj_coords;
    
    if isempty(lead_mesh_lut)
        
        lead_xyz = proj_coords{:,{'x','y','z'}};
        pial_xyz = pial_surf.vertices;
        nearest_mesh_ind = lead_to_pial(lead_xyz, [], pial_xyz, contact_radius, 0);
        n = length(nearest_mesh_ind);
        
        % Get ROIs within ROI radius of every mesh vertex mapped to this electrode
        do_lookup = @(V) roi_lut{roi_lut.d <= roi_radius & ismember(roi_lut.vertex, V), 'ROIC_mesh_ndx'};
        
        ROIC_mesh_ndxs = cellfun(do_lookup, nearest_mesh_ind, 'uniformOutput',0);
        
        
        lead_roi_lut.nearest_mesh_ind = cellfun(@(x) unique(x,'stable'), nearest_mesh_ind, 'uniformOutput',0);
        lead_roi_lut.ROIC_mesh_ndx = cellfun(@(x) unique(x,'stable'), ROIC_mesh_ndxs, 'uniformOutput',0);
        
        % add dist
        ROIC_d = [];
        for i=1:n
            dd = [];
            ROIC_mesh_ndx       = lead_roi_lut(i,:).ROIC_mesh_ndx{1};
            nearest_mesh_ind    = unique(lead_roi_lut(i,:).nearest_mesh_ind{1});
            for j = 1:numel(ROIC_mesh_ndx)
                mask = (roi_lut.ROIC_mesh_ndx == ROIC_mesh_ndx(j)) & (ismember(roi_lut.vertex, nearest_mesh_ind));
                dd(j) = min(roi_lut(mask,:).d);
            end
            ROIC_d{i} = dd(:);
        end
        lead_roi_lut.ROIC_d = ROIC_d(:);
    
    else
        lead_roi_lut = join(lead_roi_lut, lead_mesh_lut);
        
        % Filter by mesh vertices that are ROICs and within given radius
        for i = 1:height(lead_roi_lut)
            Vrow = lead_roi_lut.geodisc_mesh_ind{i};
            drow = lead_roi_lut.geodisc_dist{i};
            mask = (drow <= roi_radius) & ismember(Vrow, roi_lut.ROIC_mesh_ndx);
            ROIC_mesh_ndx{i} = Vrow(mask);
            ROIC_d{i}     = drow(mask);
        end
        
        lead_roi_lut.ROIC_mesh_ndx = ROIC_mesh_ndx(:);
        lead_roi_lut.ROIC_d = ROIC_d(:);
    end
    
    remove_columns = {'geodisc_mesh_ind','geodisc_dist'};
    keep_cols = setdiff(lead_roi_lut.Properties.VariableNames, remove_columns, 'stable');
    lead_roi_lut = lead_roi_lut(:, keep_cols);
    
    
    
end



    
    
    
    
    
    
    
    
    
    
    
    

