function util_file_transfer_horizontal(varargin)
%%
% Desc:
%   horizontal copy transfer of files
% Inputs:
% Outputs:

%%
% determine defaults
default_cell = {'from', 'test', 'parameter';
    'to','server','parameter';
    'params','fillme','parameter';
    'reset',true, 'parameter';
    'which',{},'parameter'};
inputs = util_extract_inputs(default_cell, varargin);
%%
% yo gimme some params, foo
if strcmp('params','fillme')
    error('where my params, foo?\n');
end
%%
% determine paths 
util_set_data_env(inputs.from);
path_list_from = util_flood_list(inputs.params, util_list_paths);

util_set_data_env(inputs.to);
path_list_to = util_flood_list(inputs.params, util_list_paths);

% try to copy
for sandwich = inputs.which
    try
        [try_dir try_name try_ext] = fileparts(path_list_to.(sandwich{1}));
        if ~exist(try_dir)
            mkdir(try_dir)
        end
        if findstr(path_list_to.(sandwich{1}),'[whichHemi]')
            for my_hemi = {'lh','rh'}
                xfr_stuff.from = path_list_from.(sandwich{1});
                xfr_stuff.to = path_list_to.(sandwich{1});
                xfr_stuff = util_flood_list({'whichHemi', my_hemi{1}},xfr_stuff);
                copyfile(xfr_stuff.from, xfr_stuff.to)
            end
        else
            copyfile(path_list_from.(sandwich{1}),path_list_to.(sandwich{1}));
        end
    catch
        exist(path_list_from.(sandwich{1}))
        fprintf('failed for %s\n', sandwich{1});
    end
end

% reset to default data directory after finished
if inputs.reset
    fprintf('resetting the paths via startup.m, by default\n');
    startup;
end
return