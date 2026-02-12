function [g2_idx, g1_list] = util_group2group(g1_idx, g2_list)
[g1_ind g2_ind] = util_group2ind(g1_idx,g2_list);
[g2_idx, g1_list] = util_ind2group(g2_ind, g1_ind);
return
