function util_migrate_dir(foc_from, foc_to, params)
% Desc: migrates data from focus1 to focus2
% Inputs: herpamerpaderp
% Outputs:
%%
params.subjId = 'NIH042';
foc_from = 'proc';
foc_to = 'best';
dir_list_from= util_flood_list(params, util_list_dirs('path_focus',foc_from));
dir_list_to = util_flood_list(params,util_list_dirs('path_focus',foc_to));

ignore_list = {'docs','tal','freesurfer_base','suma','inputs','intermediates','temp_db'};

field_list = fieldnames(dir_list_from);
for my_field = field_list'
    if ~ismember(my_field{1}, ignore_list)
        try
            if ~exist(dir_list_to.(my_field{1}))
                mkdir(dir_list_to.(my_field{1}))
            end
            copyfile(dir_list_from.(my_field{1}),dir_list_to.(my_field{1}));
        catch
            fprintf('aw %s\n',my_field{1})
        end
    end
end
%%
return