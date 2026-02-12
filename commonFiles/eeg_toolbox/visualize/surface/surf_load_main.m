function [my_surf] = surf_load_main(varargin)
    %% define inputs of interest
    varargin = util_ensure_varargin(varargin);
    default_cell = {
        'subjId','','parameter';
        'rootEEGdir','','parameter';
        'hemi', 'lh', 'parameter';
        'surfRez','std.141','parameter';
        'surfType', 'pial', 'parameter';
        'util2dir',0,'parameter';
        'dir','','parameter'
        };
    my_surf = util_extract_inputs(default_cell, varargin); % surface info

%     if ~isempty(my_surf.rootEEGdir)
%         my_surf.dir = fullfile(rootEEGdir, subj, 'tal/intermediates/images_2_fsSumaStd', subj);
%     end
    
    
    %% build directory if so desired 
    if my_surf.util2dir
        dir_list = util_list_dirs(inputs.rootEEGdir);
        my_surf.dir = util_fill_params({'subjId',my_surf.subjId}, dir_list.suma);
    end
    my_surf = rmfield(my_surf, 'util2dir');

    %% return based on hemisphere request
    switch my_surf.hemi
        case 'combined'
            my_surf = surf_load_both(my_surf);
            my_surf = surf_combine_hemis(my_surf);
        case 'separate'
            my_surf = surf_load_both(my_surf);
        case {'lh','rh','r','l'}
            my_surf = surf_load_surf(my_surf);
    end
    %%
return

function both_surf = surf_load_both(my_surf)
    both_surf = struct();
    % load each hemi separately
    for my_hemi = {'lh','rh'}
        my_surf.hemi = my_hemi{:};
        both_surf.(my_hemi{:}) = surf_load_surf(my_surf);
    end
return

function my_surf = surf_load_surf(my_surf)
    % build the path
    my_surf.name = surf_build_filename(my_surf.surfRez, my_surf.hemi, my_surf.surfType);
    my_surf.path = fullfile(my_surf.dir, my_surf.name);

    try
        if ~exist(my_surf.path, 'file')
            error('Expected file not found: %s', my_surf.path);
        else
            fprintf('Loading %s...\n', my_surf.name);
        end
    catch
        my_surf.surfRez = 'native';
        my_surf.name = surf_build_filename(my_surf.surfRez, my_surf.hemi, my_surf.surfType);
        my_surf.path = fullfile(my_surf.dir, my_surf.name);
        if ~exist(my_surf.path, 'file')
            error('Expected file not found: %s', my_surf.path);
        else
            fprintf('Loading %s...\n', my_surf.name);
        end
    end

    % read gifti in
    raw_surf = gifti(my_surf.path);
    my_surf.faces = round(double(raw_surf.faces));
    my_surf.vertices = double(raw_surf.vertices);
    my_surf.orientation = 'LPI';
    try my_surf.xfm = raw_surf.mat; catch; end
return

function surf_name = surf_build_filename(surfRez, hemi, surfType)
    switch surfRez
        case 'native' 
            surfRez = '';
        otherwise 
            surfRez = strcat(surfRez, '.');
    end
    surf_name = strcat(surfRez, hemi, '.',  surfType, '.gii');
return
