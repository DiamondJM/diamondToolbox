function [dset] = dset_load_main(varargin)
%% define inputs of interest
varargin = util_ensure_varargin(varargin);
default_cell = {
    'dset','sulc','parameter';
    'subjId','','parameter';
    'whichHemi', 'lh', 'parameter';
    'surfRez','std.141','parameter';
    'surfType', 'pial', 'parameter';
    'util2dir',false,'parameter';
    'dir','./','parameter'
    };
dset = util_extract_inputs(default_cell, varargin);

%%
if dset.util2dir
    dir_list = util_list_dirs;
    dset.dir = util_fill_params({'subjId',dset.subjId}, dir_list.suma);
end
dset = rmfield(dset, 'util2dir');
%% return based on which hemisphere request
switch dset.whichHemi
    case 'combined'
        dset = dset_load_both(dset);
        dset = dset_combine_hemis(dset);
    case 'separate'
        dset = dset_load_both(dset);
    case {'lh','rh'}
        dset = dset_load(dset);
end
%%
return
function both_dset = dset_load_both(dset)
both_dset = struct();
% load each whichHemi separately
for my_whichHemi = {'lh','rh'}
    dset.whichHemi = my_whichHemi{:};
    both_dset.(my_whichHemi{:}) = dset_load(dset);
end
return
function combined = dset_combine_hemis(both_dset)
combined = both_dset.lh; % initialize

for my_hemi = {'lh','rh'};
    both_dset.(my_hemi{:}).num_verts = length(both_dset.(my_hemi{:}).data);
end

combined.data = [both_dset.lh.data; both_dset.rh.data];
combined.lhIdx = [1 both_dset.lh.num_verts];
combined.rhIdx = [(both_dset.lh.num_verts+1) (both_dset.lh.num_verts + both_dset.rh.num_verts)];
combined.path     = strrep(both_dset.lh.path, 'lh', '*');
combined.name     = strrep(both_dset.lh.path,'lh','*');
return
function dset = dset_load(dset)
% build the path
dset.name = dset_build_filename(dset.surfRez, dset.whichHemi, dset.dset);
dset.path = fullfile(dset.dir, dset.name);

% read dset in
raw_dset = afni_niml_readsimple(dset.path);
dset.data = raw_dset.data;
return
function dset_name = dset_build_filename(surfRez, whichHemi, dset)
switch surfRez
    case 'native'
        surfRez = '';
    otherwise
        surfRez = strcat(surfRez,'.');
end
dset_name = strcat(surfRez, whichHemi, '.', dset,'.niml.dset');
return