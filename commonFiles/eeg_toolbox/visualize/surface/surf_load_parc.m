function my_parc = surf_load_parc(varargin)
    %% define inputs of interest
    varargin = util_ensure_varargin(varargin);
    default_cell = {
        'subjId','','parameter';
        'hemi', 'lh', 'parameter';
        'surfRez','std.141','parameter';
        'parcType','a2009s','parameter';
        'util2dir',false,'parameter';
        'dir','./','parameter'
        };
    my_parc = util_extract_inputs(default_cell, varargin); % parc info

    %% build directory if so desired 
    if my_parc.util2dir
        dir_list = util_list_dirs;
        my_parc.dir = util_fill_params({'subjId',my_parc.subjId}, dir_list.suma);
    end
    my_parc = rmfield(my_parc, 'util2dir');

    %% return based on hemisphere request
    switch my_parc.hemi
        case 'combined'
            my_parc = parc_load_both(my_parc);
            my_parc = parc_combine_hemis(my_parc);
        case 'separate'
            my_parc = parc_load_both(my_parc);
        case {'lh','rh'}
            my_parc = parc_load_parc(my_parc);
    end
end

function both_parc = parc_load_both(my_parc)
    both_parc = struct();
    % load each hemi separately
    for my_hemi = {'lh','rh'}
        my_parc.hemi = my_hemi{:};
        both_parc.(my_hemi{:}) = parc_load_parc(my_parc);
    end
end

function my_parc = parc_load_parc(my_parc)
    % build the path
    my_parc.name = parc_build_filename(my_parc.surfRez, my_parc.hemi, my_parc.parcType);
    my_parc.path = fullfile(my_parc.dir, my_parc.name);

    % extract info from file, read it in as a by-vertex mapping from label to vert
    [vert_idx_ind label_ref_ind lut] = afni_niml_annotation(my_parc.path); % read it in
    vert_idx_ind = vert_idx_ind + 1;
    % convert to by-label mapping from label to vert list
    [label_ref vert_idx] = util_ind2group(label_ref_ind, vert_idx_ind); % gather into vertex list
    label2vert = table(label_ref, vert_idx); 

    % reformat the lut table
    r = lut.table(:,1);
    g = lut.table(:,2);
    b = lut.table(:,3);
    label_ref = lut.table(:,5); %<-- mike commented out 06/28/17
    label_name = lut.struct_names;
    label_idx = (1:length(label_name))';
    lut_format = table(r,g,b,label_idx,label_ref,label_name);

    % and combine info into final label2vert table of interest
    label2vert = innerjoin(label2vert, lut_format);
    if height(label2vert) < 2
        keyboard
        [label_ref vert_idx] = util_ind2group(label_ref_ind, vert_idx_ind); % gather into vertex list
        label2vert = table(label_ref, vert_idx); 
        MIN_NODES = 50;
        % make a random color table?

    end

    % output!
    my_parc.numLabels = height(label2vert);
    my_parc.numVerts = max([vert_idx{:}]);
    my_parc.lut = label2vert;
end

function parc_name = parc_build_filename(surfRez, hemi, parcType)
    switch surfRez
        case 'native' 
            surfRez = '';
        otherwise 
            surfRez = strcat(surfRez, '.');
    end
    parc_name = strcat(surfRez, hemi, '.aparc.',  parcType, '.annot.niml.dset');
end

function combined = parc_combine_hemis(both_parc)
    combined = both_parc.lh; % initialize

    % figure out right hemi stuff
    temp.rh = both_parc.rh;
    temp.rh.lut.label_idx = temp.rh.lut.label_idx + height(both_parc.lh.lut);
    temp.rh.lut.vert_idx = cellfun(@(x) x + both_parc.lh.numVerts, temp.rh.lut.vert_idx, 'UniformOutput',false);

    % concatenate these tables
    combined.lut = [both_parc.lh.lut; temp.rh.lut];

    % and grab meta indexing info
    combined.lhVertIdx = [1 both_parc.lh.numVerts];
    combined.numVerts = both_parc.lh.numVerts + both_parc.rh.numVerts;
    combined.rhVertIdx = [(both_parc.lh.numVerts+1) (both_parc.lh.numVerts + both_parc.rh.numVerts)];

    combined.lhLabelIdx = [1 both_parc.lh.numLabels];
    combined.numLabels = both_parc.lh.numLabels + both_parc.rh.numLabels;
    combined.rhLabelIdx = [both_parc.lh.numLabels+1 combined.numLabels];

    combined.path     = strrep(both_parc.lh.path, 'lh', '*');
    combined.name     = strrep(both_parc.lh.path,'lh','*');
end
