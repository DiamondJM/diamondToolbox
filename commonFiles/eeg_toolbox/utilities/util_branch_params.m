function my_table = util_branch_params(field_name, opt_list, my_table)
% Desc: 
%   expand out table by field name and options list
% Inputs:
%   field_name: name of field to add to table
%   opt_list: list of options to expand out by
%   my_table: either an empty matrix ([]) or another parameter table that
%             was previously iterated over
% % Example code:
% clear; clc;
% 
% merp = {'hi';'bye';'blah'};
% my_table = table(merp)
% 
% field_name = 'wassup';
% opt_list  = {'hi','this','works'}; % to show that this works with cell arrays
% opt_list = [1 2 3]; % to show that this works with matrices
% my_table = util_branch_params(field_name, opt_list, my_table)

%%
if isempty(my_table)
    num_rows = 1;
else
    num_rows = height(my_table);    
end

% get number of opts and rows
num_opts = length(opt_list);

% repeat options for each current parameter row
opt_col = arrayfun(@(x) repmat(x,num_rows,1), opt_list, 'UniformOutput',false)';
opt_col = [opt_col{:}];
opt_col = opt_col(:);

if isempty(my_table)
    my_table = table(opt_col);
    my_table.Properties.VariableNames = {field_name};
else
    % repeat parameter rows for each option
    my_table = repmat(my_table,num_opts,1);
    
    % now merge opt_col and my_table
    my_table.(field_name) = opt_col;
end
return