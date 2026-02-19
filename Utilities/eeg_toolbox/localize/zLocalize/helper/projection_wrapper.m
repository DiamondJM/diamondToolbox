function [ndx_xyz_dural, ndx_xyz_nosnap, ndx_d_dural] = projection_wrapper(coord_ndx_xyz, anchor_ndx_xyz, duralSurfVerts, k_params, neighbor_mask, fmincon_info_file)
    % projection_wrapper
    %
    % Usage: [ndx_xyz_dural, ndx_xyz_nosnap, ndx_d_dural] = projection_wrapper(coord_ndx_xyz, anchor_ndx_xyz, duralSurfVerts, configVec, varargin)
    %
    % Inputs:
    %   coord_ndx_xyz     - n x 4 (id, x, y, z) CT-registered points
    %   anchor_ndx_xyz  - n x 4 (id, x, y, z) Anchor points
    %   duralSurfVerts  - m x 3 (x,y,z) dural surface points
    %   neighbor_mask   - binary adjacency matrix squareform('vector') form to define neighbors of coord_ndx_xyz
    %   k_params        - structure of parameters with fields:
    %                       e_disp_weight
    %                       e_fit_weight
    %                       e_anch_weight
    %                       e_def_weight
    %
    % Optional Input:
    %   fmincon_info_file - filepath/name of a log file to store fmincon log
    %
    % Outputs:
    %   ndx_xyz_dural   - n x 4 (id, x, y, z) locations after optimization, snapped to nearest dural point
    %   ndx_xyz_nosnap  - n x 4 (id, x, y, z) locations after optimization without snapping
    %   ndx_d_dural     - n x 4 (id, x, y, z) distance of of final (nosnap) locations to nearest dural point

    %  Copyright (C) 2017 Mike Trotta
    
    CLIP_DURA_PERCENT = 3; % for efficiency, limit "nearest dural surface" to X percent of full dural surface

    % parse input
    duralSurfVerts = double(duralSurfVerts);
    if nargin < 6, fmincon_info_file = []; end
    
    % originally registered coordinates
    xyz_reg = coord_ndx_xyz(:,2:4);
    xyz_reg_copy = xyz_reg;

    % snap the anchors to the dural surface
    if isempty(anchor_ndx_xyz)
        anchor_ndx = [];
    else
        % gets nearest dural point to each anchor
        anchor_ndx      = anchor_ndx_xyz(:,1);
        anch_xyz        = anchor_ndx_xyz(:,2:4);
        d_anch2dural    = pdist2(anch_xyz, duralSurfVerts);
        [~, minndx]     = sort(d_anch2dural, 2);
        dural_anch_xyz  = duralSurfVerts(minndx(:,1),:);

        % now overwrite anchor point values with the nearest dural point
        dural_anchor_ndx_xyz = [anchor_ndx_xyz(:,1) dural_anch_xyz];
        anchor_ndx_xyz = dural_anchor_ndx_xyz;
    end
    
    % set projection constants
    num_pts = length(coord_ndx_xyz);
    proj_params.e_disp_weight            = k_params.e_disp_weight;
    proj_params.e_fit_weight             = k_params.e_fit_weight * ones(num_pts,1);
    proj_params.e_fit_weight(anchor_ndx) = k_params.e_anch_weight;
    proj_params.e_def_weight             = k_params.e_def_weight;

    %
    proj_params.neighborMask = neighbor_mask;
    
    % originally registered geometry as distance matrix
    proj_params.dvec_old = pdist(xyz_reg);

    % Gets the closest X percent of dural points to each registered point
    close_dural_pnts = get_nearest_pts(xyz_reg, duralSurfVerts, CLIP_DURA_PERCENT/100);

    % define objective function
    objfun = @(xyz) eTotalFun(xyz, xyz_reg, proj_params, anchor_ndx_xyz, close_dural_pnts);

    % set fmincon options
    maxFunEvals = 50 * 1000;
    maxIter = 150;      
    tolFun = 0.05;
    options = optimset('Algorithm','interior-point','Display','iter-detailed','MaxFunEvals', maxFunEvals,'MaxIter',maxIter,'TolFun',tolFun,'UseParallel',1);
    % display iter-detailed, off, final-detailed

    % run optimization
    [xyz_new,fval,exit_flag,info] = fmincon(objfun, xyz_reg_copy, [],[],[],[],[],[], [], options);
    if ~isempty(fmincon_info_file)
        save(fmincon_info_file, 'xyz_new', 'fval', 'exit_flag', 'info', 'proj_params', 'maxFunEvals','maxIter','tolFun');
    end

    % This part sets dural_xyz_new to nearest dural point (out of original dural subset) to final xyz locations
    d_new2dural = pdist2(xyz_new, duralSurfVerts);          % dist of xyz to every dural pt
    [sortd_new2dural, new_minndx] = sort(d_new2dural, 2);    % d_new2dural(k,new_minndx(k,1)) = sortd_new2dural(k,1); k'th electrode's shortest distance to dural
    dural_xyz_new = duralSurfVerts(new_minndx(:,1),:);      % new_minndx(k,1) gives index of k'th electrode's nearest dural point

    % output
    ndx_xyz_dural  = [coord_ndx_xyz(:,1) dural_xyz_new];
    ndx_xyz_nosnap = [coord_ndx_xyz(:,1) xyz_new];
    ndx_d_dural    = [coord_ndx_xyz(:,1) sortd_new2dural(:,1)];

end

function pts = get_nearest_pts(reg_xyz, duralVerts, percPts)
    % Finds the closest X percent of dural surface points closest to each
    % electrode. For each reg_xyz point, returns a set of those closest points

    numPts = round(percPts*length(duralVerts(:,1))); % number of points

    % for each electrode coordinate, select the set of all possible points
    % to search among for the closest face in the findSurface call of fmincon
    nxyz = length(reg_xyz);
    pts = cell(nxyz, 1);
    for i = 1 : nxyz

        % find the closest X dural points
        pt_i = reg_xyz(i, :);
        d = pdist2(pt_i, duralVerts, 'euclidean');
        [~,ndx] = sort(d); 

        pts{i} = duralVerts(ndx(1:numPts),:);
    end
end

function eTotal = eTotalFun(xyz, xyz_old, proj_params, anchor_idx_xyz, close_dural_pnts, debug)
    % Calculates energy functions: E_deformation, E_displacement, E_fit
    % close_dural_pnts is closest envelope points to each electrode
    % e_coda is E_total

    if ~exist('debug', 'var')
       debug = 0;
    end

    % load paramaters of interest
    dvec_old = proj_params.dvec_old;       % registered pairwise distances
    Nmask = logical(proj_params.neighborMask); % neighbor mask
    k_disp = proj_params.e_disp_weight;    % ct registered
    k_fit = proj_params.e_fit_weight;      % projection/anchor weight vector
    k_def = proj_params.e_def_weight;      % deformation / geometry

    % get e_deformation term
    dvec = pdist(xyz);

    % apply neighbor mask
    nbr_dvec_old = dvec_old(Nmask);
    nbr_dvec = dvec(Nmask);

    % calculate mean of squared distance difference
    e_def_vec = (nbr_dvec - nbr_dvec_old) .^ 2;
    e_def =  mean(k_def * e_def_vec);

    % get e_displacement term
    e_disp_vec = sum( (xyz - xyz_old) .^ 2, 2 ); % ||x - x'||^2 for each difference
    e_disp = mean(k_disp * e_disp_vec);

    % get e_fit (anchor and projection)
    xyz_projected = findSurface(xyz, close_dural_pnts);
    if debug
        assignin('base','xyz_proj',xyz_projected)
    end

    % overwrite anchor electrodes with their given anchor coordinates
    if ~isempty(anchor_idx_xyz)
        xyz_projected(anchor_idx_xyz(:,1),:) = anchor_idx_xyz(:,2:4);
    end

    e_fit_vec = sum( (xyz - xyz_projected) .^ 2, 2 ); % ||x - x'||^2 for each difference
    e_fit = mean(k_fit .* e_fit_vec);

    % add to get e_total
    eTotal = e_def + e_disp + e_fit;

    if debug
        k_fit = mode(k_fit); % probably less anchors than points
        k_anch = setdiff(proj_params.e_fit_weight, k_fit); 
        if isempty(k_anch), k_anch = k_fit; end
        if isempty(k_anch), k_anch = 0; end
        fprintf('%d =\t%d\t%d\t%d \t(disp %d, fit %d + anch %d, def %d)\n', eTotal, e_disp,e_fit,e_def, k_disp,k_fit,k_anch,k_def);
    end

end


function proj_pts =  findSurface(xyz, all_dural_pts)
	% Finds the closest point on the surface. 
    % Given our original set of points closest to orignally registered points,
    % find the closest dural point to the current xyz point

	proj_pts = zeros(length(xyz), 3);
    
    % get the nearest dural point for each xyz
	for k = 1:size(xyz,1)
	    close_durals = all_dural_pts{k};
	    x_temp = xyz(k,:);
	    dist = pdist2(x_temp,close_durals);
	    [~,idx] = sort(dist,'ascend');
	    proj_pts(k,:) = close_durals(idx(1),:);
	end
end
