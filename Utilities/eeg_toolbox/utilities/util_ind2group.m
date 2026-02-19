function [uni_group, form_idx] = util_ind2group(group, idx, dummy)
% nargin = 2;
% group = {'A','A','B','B','C','C','D','D'}'
% % group = [1 1 2 2 3 3 4 4]';
% idx = [1 2 3 4 5 6 7 8]';
%%
if nargin == 2
    transposer = 1;
else 
    transposer = 0;
end
%%
[group_idx uni_group] = findgroups(group);
if transposer
    form_idx= splitapply(@(x) {x'}, idx, group_idx);
%     form_idx = splitapply(@(x) length(x),idx,group_idx);
else
    form_idx= splitapply(@(x) {x}, idx, group_idx);
end
%%
% uni_group = unique(group);
% form_idx = {};
% 
% if length(uni_group) > 1000
%     warning('you have a crap ton of groups, you sure this is right?\n');
% end
% if length(uni_group) > 10000
%     error('more than 10k elements, im exiting this shiz');
%     uni_group = 'FAIL';
%     form_idx = 'FAIL';
% else
%     for i_group = 1:length(uni_group)
%         my_group = uni_group(i_group);
%         if transposer
%             form_idx{i_group,1} = idx(ismember(group,my_group),:)';
%         else
%             form_idx{i_group,1} = idx(ismember(group,my_group),:);
%         end
%     end
% end
%%
return
