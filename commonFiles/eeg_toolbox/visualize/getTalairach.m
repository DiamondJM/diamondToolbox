function [tal_labels, dist, fsaverage] = getTalairach(vertex, varargin)
% GETTALAIRACH returns the Talairach atlas of a given std141 "HD" vertex
%
% INPUTS:
%   vertex - node (1 to 198,812) from the std 141 tesselation
%
% OPTIONAL:
%   within=X - Filter atlases returned to within X (default=2) mm
%   fsaverage=surf - pass in the fsaverage surface to avoid loading it here
%
% OUTPUTS:
%   tal_labels - cell array of labels
%   within - distance of atlas to node
%   fsaverage - fsaverage surface (left hemisphere)
%
% DESCRIPTION
%   This is a rough approximation implemented as follows:
%    1) Obtain xyz-coordinate of given vertex on FSAverage brain (MNI
%       coord)
%    2) Use mni2tal *rough* transform to talairach space
%    3) Parse afni whereami command output of xyz for atlas labels
%   
%
% REVISION HISTORY
%   11/16 MST - Created
% 

ip = inputParser;
ip.addParameter('within', 2);
ip.addParameter('fsaverage', []);
ip.parse(varargin{:});
within = ip.Results.within;
fsaverage = ip.Results.fsaverage;

tal_labels = {};
dist = [];

% load fsaverage surface
if isempty(fsaverage)
    bd = braindata;
    bd.loadFSAverage();
    fsaverage = bd.surf.lh;
else
    if ~ismember('vertices', fieldnames(fsaverage))
        fsaverage = fsaverage.lh;
    end
end

xyz_mni = fsaverage.vertices(vertex,:);
xyz_tal = mni2tal(xyz_mni);


atlas_table = afni_extract_whereami('xyz', xyz_tal, 'atlas', 'TT_Daemon', 'orientation', 'lpi', 'run',true);

% fail
if isempty(atlas_table), return; end

atlas_table = atlas_table(atlas_table.within_mm <= within, :);
tal_labels = atlas_table.label;
dist = atlas_table.within_mm;

end