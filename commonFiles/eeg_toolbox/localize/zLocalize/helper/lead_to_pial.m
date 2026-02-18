function pial_nn = lead_to_pial(leads_xyz, dural_xyz, pial_xyz, r, snap_to_dura_first)
    % Desc:
    %   attribute nodes on pial to dural defined lead locations
    %
    %
    %   Given a XYZ locations of target points (leads_xyz), for each electrode, find
    %   the nearest point in dural_xyz, and from this point find a set of all other
    %   dural_xyz points within euclidean distance r. Now for each of these
    %   sets of points, find each point's nearest point in pial_xyz
    %
    %   If snap_to_dura_first is true, then just find the nearest pial points
    %   to a spherical-modeled electrode 
    %
    %
    % Inputs:
    %   leads_xyz: leads locations, assumed to be within (or near) dural set
    %   dural_xyz: dural vertices (required if snap_to_dura_first)
    %   pial_xyz:  pial vertices
    %   r:         radius of dural growth
    %   
    % Optional Inputs:
    %   snap_to_dura_first: If false, radius centered at given leads_xyz. If
    %                       true (default), radius centered at point in
    %                       dural_xyz nearest to given leads_xyz.
    %
    % Outputs:
    %   pial_nn: pial nearest neighbors as cells of arrays of indices into pial_xyz surface (1 cell per leads_xyz row)
    %
    % REVISION HISTORY
    %   03/17 MST - Change default behavior to NOT snap to dura initially
    % 
    if nargin < 5
        snap_to_dura_first = false;
    end

    nleads = length(leads_xyz(:,1));
    pial_nn = cell(nleads, 1);

    % get closest dural indices within a sphere

    % snap to dura
    if snap_to_dura_first
        leads_xyz = dural_xyz(knnsearch(dural_xyz, leads_xyz),:); 
        dmat = pdist2(leads_xyz, dural_xyz);
        filter_dmat = dmat <= r;

        dural_nodes = cell(nleads, 1);
        for ilead = 1:nleads
            dural_nodes{ilead} = find(filter_dmat(ilead,:));
        end

        % map dural nodes to nearest neighbor pial nodes
        for ilead = 1:nleads
            pial_nodes = [];
            my_dural_xyz = dural_xyz(dural_nodes{ilead}, :);
            pial_nodes = knnsearch(pial_xyz, my_dural_xyz);
            pial_nn{ilead} = unique(pial_nodes);     
        end
        
    else
        
        % model each electrode as a sphere point-cloud
        npoints = 21;
        [x,y,z] = sphere(npoints-1);
        baseSphere = permute([x(:),y(:),z(:)], [3 2 1]);
        spheres = repmat(baseSphere, nleads, 1); % spheres will be lead X [x,y,z] X npoints^2
        
        if contains(lower(which('knnsearch')), 'ielvis')
            error('Remove bad knnsearch from path: ielvis');
        end
        
        parfor i = 1:nleads
            mySphere = spheres(i,:,:) .* r + leads_xyz(i,:);
            spheres(i,:,:) = mySphere;
            contact_xyz = squeeze(mySphere)'; % npoints x 3
            pial_nn{i} = unique(knnsearch(pial_xyz, contact_xyz));
        end
        
        
        
    end

    N_PIAL_WARN_NODE = 100;
    if any(cellfun(@numel, pial_nn) > N_PIAL_WARN_NODE)
        warning('lead_to_pial: large numbers of pial nodes found for at least 1 lead (bad which(knnsearch)?) ');
    end

end
