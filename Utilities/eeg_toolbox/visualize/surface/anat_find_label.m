function vert_match = anat_find_label(match_struct, targ_name)
% Desc: for a given targ_name, find vertex indices that match
% Inputs:
%   match_struct.vert_table
%       vert_idx
%       label_idx
%   match_struct.label_table
%       r,g,b
%       label_idx
%       label_name
% Outputs:
%   vert_match: vert_idx that match 
% get vertices that correspond to given target

% targ_name --> label_idx --> vert_idx
rel_idx = strcmp(match_struct.label_table.label_name,targ_name);
abs_idx = match_struct.label_table.label_idx(rel_idx);
abs_match = match_struct.vert_table.label_idx == abs_idx;
vert_match = match_struct.vert_table.vert_idx(abs_match);
return

