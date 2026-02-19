function g_table = util_table_gather(my_table, varargin)
%%
default_cell = {'name','val','parameter'; 
                'by','na','parameter'}; % we need this to be required
params = util_extract_inputs(default_cell, varargin);
%%
by_table     = my_table(:,params.by);
not_by_table = my_table; not_by_table.(params.by) = [];
%%
for i_row = 1:height(my_table)
    by_data = by_table{i_row,params.by}{:};
    g_data = repmat(not_by_table(i_row,:), length(by_data),1);
    g_data.(params.name) = by_data;
    if i_row == 1
        g_table = g_data;
    else
        g_table = [g_table; g_data];
    end
end
%%
return



% this function will operate as follows
%   for each row of the table
%       repeat non by variables
%       and split out the by variables


% so quite literally...
% all i have to do is say...

