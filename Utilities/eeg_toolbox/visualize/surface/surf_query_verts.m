function [vert_match] = surf_query_verts(label, vert_arr, targ)
% Desc:
% Inputs:
%   label: array of labels
%   vert_arr: vertex array associated with labels
%   targ: target name or string
% Outputs:
%   vert_match: matching vertex
%%
% label = {'A','B','C'}';
% targ = 'A';
% vert_arr = {[1 2]; [3 4]; [5 6]};
vert_match = vert_arr(ismember(label, targ)); % there should only be one
if iscell(vert_arr)
    vert_match = vert_match{:};
end
%%
return