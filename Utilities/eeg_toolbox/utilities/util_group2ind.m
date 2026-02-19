
function [form_group, form_idx] =  util_group2ind(group, idx)
% group = {'A';'B';'C'};
% idx = {[1 2 3]; [4 5 6]; [7 8 9]}
% form_group = group([]);
% form_idx = [];
% for i_group = 1:length(group)
%   my_group = group(i_group);
%   my_list = idx{i_group};
%   for my_idx = my_list
%     form_group = [form_group; my_group];
%     form_idx = [form_idx; my_idx];
%   end
% end

%%
% for each group
%   repmat group for length of each of these
%   and then concatenate all that shiz

idx_length = cellfun(@length, idx, 'UniformOutput',true);
group = group(:)'; group_list = cell(length(group),1);
for i_group = 1:length(group)
    group_list{i_group} = repmat(group(i_group), 1, idx_length(i_group));
end
form_group = [group_list{:}]';
%%
form_idx = [idx{:}]';
%%
return
