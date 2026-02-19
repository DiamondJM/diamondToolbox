function [lead_mesh_lut] = lead_to_mesh_grow(pial_surf, contact_r, r, t_lead_xyz, hemi, filename_save)
    % Takes a set of channels and grows a geodesic disc around them on the pial surface
    % 
    %
    %
    % The purpose of this is to create a file which stores a region (and distance) around each electrode
    % This information can then be used to
    %   1) create subject-specific region plots
    %   2) feed into grow_ROIs to reduce the set of ROICs to only those relevant (ie near an electrode)
    %
    % Inputs
    %   pial_surf   - single hemisphere surface as struct with .vertices, .faces fields
    %   contact_r   - radius of an electrode
    %   r           - radius of growth (geodesic disc)
    %   t_lead_xyz  - table of *this* hemisphere's leads with columns: chanName, x, y, z
    %   hemi        - 'lh' or 'rh'
    %   filename_save [OPTIONAL] - if given, create/update the mat file with the given info
    %
    % Output
    %   lead_mesh_lut with fields FILTERED BY HEMISPHERE
    %       chanName
    %       whichHemi
    %       nearest_mesh_ind
    %       geodisc_mesh_ind
    %       geodisc_dist
    %
    %   leads_mesh_lut may be saved to file given by filename
    %
    % Revision History
    %   01/18 MST   - Created
    %   02/18 MST   - Made hemi optional and t_lead_xyz accept just the xyz (for cylinder compatibility)
    
    if nargin < 6 || isempty(filename_save)
        filename_save = [];
    end
    if nargin < 5 || isempty(hemi)
        hemi = '';
    end
    

    % Map electrode to pial nodes (via sphere model)
    pial_surf.vertices = single(pial_surf.vertices);
    pial_xyz = pial_surf.vertices;
    use_dura = false;
    dura     = [];
    if istable(t_lead_xyz)
        lead_xyz = t_lead_xyz{:,{'x','y','z'}};
    else
        lead_xyz = t_lead_xyz;
    end
    
    nearest_mesh_ind = lead_to_pial(lead_xyz, dura, pial_xyz, contact_r, use_dura);
    n  = length(nearest_mesh_ind);
    
    
    % Add lead and its nearest mesh points
    lead_mesh_lut = struct();
    lead_mesh_lut.whichHemi = repmat({hemi}, n, 1);
    lead_mesh_lut.nearest_mesh_ind = nearest_mesh_ind;
    if istable(t_lead_xyz)
        lead_mesh_lut.chanName = t_lead_xyz.chanName;
    else
        lead_mesh_lut.chanName = repmat({''}, n, 1);
    end
    
    % Find and add a the geodesic disc around each lead
    for i = 1 : n
        V_center = nearest_mesh_ind{i};
        tgrow = mesh_grow_centers(pial_surf, V_center, r); % note t_grow has columns "vertex" and "d"
        [~,i_dsort] = sort(tgrow.d);
        tgrow = tgrow(i_dsort,:);
        [~,i_vuniq] = unique(tgrow.vertex, 'stable'); % Take the nearest one if there are duplicates
        tgrow = tgrow(i_vuniq,:);
        
        geodisc_mesh_ind{i} = double(tgrow.vertex(:));
        geodisc_dist{i} = tgrow.d(:);
    end
    lead_mesh_lut.geodisc_mesh_ind  = geodisc_mesh_ind(:);
    lead_mesh_lut.geodisc_dist      = geodisc_dist(:);
    lead_mesh_lut = struct2table(lead_mesh_lut);
    
    % Update if previous
    if ~isempty(filename_save)
        if exist(filename_save, 'file')
            data = load(filename_save);
            lead_mesh_lut_orig = data.lead_mesh_lut;
            addFromOld = ~ismember(lead_mesh_lut_orig.chanName, lead_mesh_lut.chanName);
            
            % add any old entries which aren't in the new one to the new one:
            lead_mesh_lut = cat(1, lead_mesh_lut_orig(addFromOld,:), lead_mesh_lut);
        end
        save(filename_save, 'lead_mesh_lut');
    end

    lead_mesh_lut = lead_mesh_lut(strcmpi(lead_mesh_lut.whichHemi, hemi), :);
end