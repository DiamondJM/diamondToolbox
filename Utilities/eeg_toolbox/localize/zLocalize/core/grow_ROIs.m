function roi_mesh_d_lut = grow_ROIs(subj_surf, std_ROICs, radius, varargin)
% grow_ROIs 
%
% Usage: roi_mesh_d_lut = grow_ROIs(subj_surf, std_ROICs, radius, ...)
%
% Inputs:
%   subj_surf       - Pial surface. This is a struct which must have the following fields:
%       vertices      - n x 3 array of xyz-coordinates where n is the number of vertices in the standard SUMA mesh
%       faces         - m x 3 array of mesh faces with values in [1,n]
%       
%   std_ROICs       - Array of ROI centers stored as indices into the surface mesh vertices (ie in [1,n])
%   radius          - Radius to grow geodesic ROI discs (millimeters)
%   
% Input key-values:
%   log_file        - path to log file. If not supplied, print output to console.
%   
% Output: 
%   roi_mesh_d_lut  - Look-up table storing distances from each ROI center to mesh indices within a radius with these columns:
%       vertex        - Index into pial surface mesh (standard SUMA-resampled)
%       d             - Distance from ROIC to vertex (mm)
%       ROIC_mesh_idx - Index into surface mesh of the ROIC (in [1,n])
%       ROIC_idx      - Index 1 to k of this row's ROI center among all k ROI centers
%
% File Output:
%   
%
% Description:
%   Places the standard ROICs from the standard brain onto the subject brain, then grows ROIs
%   geodesically from the centers for a given radius and stores mesh vertices' distances to the ROICs
%   

% Revision History:
%   03/17 - MST

%  Copyright (C) 2017 Mike Trotta

    
    % Input
    ip = inputParser;
    ip.addParameter('log_file', []);
    ip.parse(varargin{:});
    log_file = ip.Results.log_file;
    
    required_fields = {'vertices', 'faces'};
    assert(length(fieldnames(subj_surf)) >= length(intersect(fieldnames(subj_surf), required_fields)), ...
        'subj_surf must have fields: vertices, faces');

    if ~isempty(log_file) && ~exist(log_file, 'file'), error('Log file does not exist: %s', log_file); end
        
    roi_mesh_d_lut = mesh_grow_centers(subj_surf, std_ROICs, radius, log_file);
    
end
