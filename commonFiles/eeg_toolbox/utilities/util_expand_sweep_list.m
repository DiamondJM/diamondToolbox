function params_table  = util_expand_sweep_list(sweep_list, params_table)
field_list = fieldnames(sweep_list);
num_fields = length(field_list);
for i_field = 1:num_fields
    my_field = field_list(i_field); my_field = my_field{:};
    opt_list = sweep_list.(my_field);
    params_table = util_branch_params(my_field, opt_list, params_table);
end
return