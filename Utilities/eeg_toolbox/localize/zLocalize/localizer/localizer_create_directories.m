function dirs = localizer_create_directories(subj, rootEEGdir, numCTs)
% dirs = localizer_create_directories(subj, rootEEGdir, numCTs)
%
%   Inputs
%       subj       - code name of subject
%       rootEEGdir - rootEEGdir should be folder containing subject folders
%       numCTs     - Optional, default 1. Useful in the case of a rearrangement of electrodes for a single subject
%
% This file is intended to be a template for your own localization pipeline. It serves
% to create and return and structure of directories under the root data directory.
%
%
% History
%   01/18 MST   - Moved anchors, atlas, and roi folders up to tal/ directory
%

%  Copyright (C) 2017 Mike Trotta

    if nargin < 3, numCTs = 1; end
    assert(numCTs == 1 || numCTs == 2, 'Localization pipeline not built for more than 2 CT files');
    
    
    % create subject directory if it doesn't exist
    tal  = fullfile(rootEEGdir, subj, 'tal');
    zloc = fullfile(rootEEGdir, subj, 'tal/zloc');
    if ~exist(zloc, 'dir'), mkdir(zloc); end

    
    
    dirs.ct_1       = makedir(fullfile(zloc,'CT_1'));
    dirs.ct_1_dicom = makedir(fullfile(zloc,'CT_1','dicoms'));
    dirs.ct_1_xfm   = makedir(fullfile(zloc,'CT_1','transform'));
    dirs.ct_2       = makedir(fullfile(zloc,'CT_2'));
    dirs.ct_2_dicom = makedir(fullfile(zloc,'CT_2','dicoms'));
    dirs.ct_2_xfm   = makedir(fullfile(zloc,'CT_2','transform'));
    dirs.fs         = makedir(fullfile(zloc,'freesurfer'));
    dirs.fs_subj    = makedir(fullfile(zloc,'freesurfer',subj));
    dirs.mr_pre     = makedir(fullfile(zloc,'mr_pre'));
    dirs.mr_t2      = makedir(fullfile(zloc,'mr_pre_t2'));
    dirs.mr_pre_dicom = makedir(fullfile(zloc,'mr_pre','dicoms'));
    dirs.mr_t2_dicom = makedir(fullfile(zloc,'mr_pre_t2','dicoms'));
    dirs.mr_resect  = makedir(fullfile(zloc, 'mr_resect'));
    dirs.mr_implant = makedir(fullfile(zloc, 'mr_implant'));
    dirs.working    = makedir(fullfile(zloc,'working'));
    dirs.metrics    = makedir(fullfile(zloc, 'metrics'));
    dirs.anchors    = makedir(fullfile(zloc,'anchors'));
    dirs.alice      = fullfile(zloc,'ALICE');
    dirs.logs       = makedir(fullfile(zloc,'logs'));
    

    % tal/ roi,anchors,atlas folders
    if exist(fullfile(zloc, 'roi'), 'file')
        system(sprintf('mv %s %s', fullfile(zloc,'roi'), tal));
    end
    if exist(fullfile(zloc, 'atlas'), 'file')
        system(sprintf('mv %s %s', fullfile(zloc,'atlas'), tal));
    end
    
    dirs.roi        = makedir(fullfile(tal,'roi'));
    
    dirs.atlas      = makedir(fullfile(tal,'atlas'));
        
        
end

function dirname = makedir(dirname)
    if ~exist(dirname, 'dir')
        mkdir(dirname);
    end
end