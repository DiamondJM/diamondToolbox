function out_list = util_floods_list(params_cell, my_list)
%   Desc:
%       flood a list with parameters of interest
%   Inputs: 
%       params_cell: name-value parameter cell array, in tabular format
%   Outputs:
%       my_list: list with parameters replaced 

%%
switch class(my_list)
    case 'struct'
        myfun = @(li) structfun(@(my_field) util_fill_params(params_cell, my_field), li, 'UniformOutput',0);
% WORK IN PROGRESS BELOW::::::
        %     case 'cell'
%         myfun = @(li) cellfun(@(my_field) util_floods_list(params_cell, li), 'UniformOutput', 0);
%     case 'string'
%         myfun = @(li) util_fill_params(params_cell, {li});
% I WANT THIS TO WORK :( - john
end
out_list = myfun(my_list);
return