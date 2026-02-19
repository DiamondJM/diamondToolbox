function [proj_coords, dural_coords, dural_dist] = projection(mr_coords, envelope_surf, anchor_coords, element_info, varargin)
% projection 
%
% Usage: proj_coords = PROJECTION(mr_coords, envelope_surf, anchor_coords, element_info
%
% Inputs:
%   mr_coords       - table of xyz-coordinates relative to pre-operative MRI
%   envelope_surf   - A struct with faces and vertices fields representing the cortical envelope mesh
%   anchor_coords   - table of anchor xyz-coordinates for subset of electrodes in pre-operative MRI space.
%   element_info    - A table which describes the layout / geometry of strip/grid hardware
%                     Each row is a hardware piece. And columns contain at least the following fields
%                     - bigDim (numeric, use 1 for strips)
%                     - smallDim (numeric)
%                     - isCut (string of matlab array e.g. '[]' or '[1 2 9 10]'). Will be used with an eval().
%                       Use this for modified grids. If contacts 8,16,24,32 are cut, value would be '[8 16 24 32]'
% Optional Input:
%
% Input key-values:
%   log_file        - path to log file. 
%   fmincon_file    - path to a .mat file which will contain fmincon data.
%   
% Output: 
%   proj_coords     - table of projected coordinate locations. This table will be of the same form as
%                     the mr_coords table
%   dural_coords    - same as proj_coords, but coords have been snapped to nearest point on dural
%   dural_dist      - table in the same form as mr_coords with a dist2dural column giving distance in mm to nearest point on dural
%
% File Output:
%        
%
% Description:
%   Given a set of XYZ coordinates of electrodes, a  cortical envelope surface, 
%   and coordinates of anchor points, solves an spring-model to project electrodes to the 
%   cortical envelope. Use this for one hemisphere at a time.
%   

% Revision History:
%   03/17 - MST
                                  
%  Copyright (C) 2017 Mike Trotta
    
    % -------------
    %%--- Setup ---
    % -------------
    
    k_params.e_disp_weight  = 1;    % 1    - distance to original CT coordinates
    k_params.e_fit_weight   = 25;   % 25   - distance to dural surface
    k_params.e_anch_weight  = 200;  % 200  - distance to anchor coordinate
    k_params.e_def_weight   = 1000; % 1000 - change in inter-electrode distances from CT coordinates
    
    proj_coords = table();
    dural_coords = table();
    dural_dist = table();
    
    % Input
    ip = inputParser;
    ip.addOptional('working_dir', pwd, @exist);
    ip.addParameter('log_file', @ischar);
    ip.addParameter('fmincon_file', @ischar);
    ip.parse(varargin{:});
    log_file = ip.Results.log_file;
    fmincon_file = ip.Results.fmincon_file;
    
    % Check input
    required_cols = {'chanName', 'x','y','z'};
    assert(istable(mr_coords) && length(required_cols) == length(intersect(required_cols, mr_coords.Properties.VariableNames)),...
        'mr_coords must be a table with chanName, x, y, and z columns');
    if ~isempty(anchor_coords)
        assert(istable(anchor_coords) && length(required_cols) == length(intersect(required_cols, anchor_coords.Properties.VariableNames)),...
            'anchor_coords must be a table with chanName, x, y, and z columns');
    end

    if isempty(mr_coords)
        warning('No mr_coords to project');
        return;
    end
    
    if ischar(envelope_surf)
        assert(exist(duralSurfVerts, 'file') > 0, 'envelope_surf must be a valid gifti filename or struct with faces and vertices fields');
        envelope_surf = gifti(duralSurfVerts);
    end
    required_fields = {'faces','vertices'};
    if sum(ismember(required_fields, fieldnames(envelope_surf))) < 2
        error('Pial and evelope surfaces must be structure with .faces and .vertices fields');
    end
    duralSurfVerts = envelope_surf.vertices;
    
    % tables to mats
    if isempty(anchor_coords)
        anchor_ndx_xyz = [];
        coords_ndx_xyz  = cat(2, [1:height(mr_coords)]', mr_coords{:,{'x','y','z'}});
    else
        uncommon_anchors = setdiff(anchor_coords.chanName, mr_coords.chanName);
        if ~isempty(uncommon_anchors)
            warning('The following anchors were not found in mr_coords: %s', strjoin(uncommon_anchors, ' '));
        end
        anchor_coords   = anchor_coords(~ismember(anchor_coords.chanName, uncommon_anchors), :);
        % these next few lines ensure anchor ndx matches that in mr_coords
        [~,temp] = ismember(mr_coords.chanName, anchor_coords.chanName);
        [anchor_ndx,~,sort_ndx] = find(temp);
        anchor_xyz      = anchor_coords{sort_ndx,{'x','y','z'}};
        coords_ndx_xyz  = cat(2, [1:height(mr_coords)]', mr_coords{:,{'x','y','z'}});
        anchor_ndx_xyz  = cat(2, anchor_ndx, anchor_xyz);
    end
    
    assert(~isempty(mr_coords), 'MR coords variable empty');
    neighbor_mask = getNeighborMask(element_info, mr_coords);
    
    % ----------------------
    %%--- Run Projection ---
    % ----------------------
    
    [ndx_xyz_dural, ndx_xyz_nosnap, ndx_d_dural] = projection_wrapper(coords_ndx_xyz, anchor_ndx_xyz, duralSurfVerts, k_params, neighbor_mask, fmincon_file);
    proj_coords = mr_coords;
    proj_coords{:,{'x','y','z'}} = ndx_xyz_nosnap(:,2:4); % output of algorithm (without any dural snap constraint!)
    dural_coords = mr_coords;
    dural_coords{:,{'x','y','z'}} = ndx_xyz_dural(:,2:4);
    dural_dist = mr_coords(:,'chanName');
    dural_dist{:,{'dist2dural'}} = ndx_d_dural(:,2);
end
