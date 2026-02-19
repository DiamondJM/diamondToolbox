function label_cell = surf_nodes2parcels(label_list, node_list_cell, targ_nodes_cell)
%%
% Desc: 
%   
% Inputs: 
%   label_list: list of labels
%   node_list_cell: list of nodes accompanying each label
%   targ_nodes_cell: list of nodes that we want to find, multiple list per
%   cell
% Outputs: 
%   label_cell = labels associated with each targ_nodes_cell
%%
ind = table();
[ind.label ind.node] = util_group2ind(label_list, node_list_cell); % break out into individual format
label_cell = cellfun(@(nodes) ind.label(ismember(ind.node, nodes)), targ_nodes_cell, 'UniformOutput',false);
%%
return