function localizer_run_template(subj, rootDir, hemisphere)
% localizer_run_template(subj, rootDir, hemisphere)
%
%   Inputs
%       subj       - codename of subject. This will be the title of a directory in rootDir.
%       rootDir    - directory where this subject's localizer directories will be written
%       hemisphere - 'r' for right or 'l' for left
%
%   Note: 
%       This file is intended to be a template for your own localization pipeline. Its main purpose
%       is to manage the directory structure and necessary input files, and to call the core
%       zLocalize files in sequence. You should tailor your own localizer_run file to your needs.
%
%
%   Details:
%       In this setup, rootDir is a directory which will store data for each subject. For each subject
%       a structure of various directories will be created. The pipeline will look for necessary inputs
%       in this newly created directory structure. If it does not find the needed input files, dialog
%       boxes will allow you to locate them, and the program will copy the files into their correct 
%       locations in [rootDir]/[subj]
%

%  Copyright (C) 2017 Mike Trotta

    % ---------- Settings ---------- %
    FS_BIN          = '/Applications/freesurfer/bin/';                  % path to freesurfer bin/ directory
    AFNI_BIN        = '/Users/trottams/abin';                           % path to afni's abin/ directory
    SLICER_BINARY   = '/Applications/Slicer.app/Contents/MacOS/Slicer'; % path to the 3DSlicer binary file
    N_ROIC          = 2400; % (coarse=600, fine=2400)
    R_ROI           = 9;    % (course=14, fine=9)
    R_CONTACT       = 1.5;  % electrode contact radius
    USE_SYM_LINK    = 0;    % Whether or not to use a symbolic link for freesurfer
    
    assert(exist(rootDir,'dir') > 0, 'Root directory does not exist: %s\n', rootDir);
    
    fprintf('\n\nRunning localizer_run_template, a template zLocalize pipeline.\n');
    
    fprintf('Creating a directory structure in which zLocalize will store all necessary files: %s\n', fullfile(rootDir, subj));
    fprintf('\tcreating...');
    locDirs = localizer_create_directories(subj, rootDir);
    fprintf('done!\n');    
    
    fname_default_ct_coord      = fullfile(locDirs.ct_1, 'coords.csv');
    fname_default_element_info  = fullfile(locDirs.docs, 'element_info.csv');
    fname_default_ROIC          = which(sprintf('zLocalize_ROI_centers_%sh.mat', hemisphere));
    fname_default_fs_subj       = fullfile(locDirs.fs, subj);
    fname_default_mr_nii        = fullfile(locDirs.mr_pre, 'mr_pre.nii');
    fname_default_ct_nii        = fullfile(locDirs.ct_1, 'ct_implant.nii');
    
    % locate / copy into rootDir basic files
    ct_coords       = readtable(uigetfile_if_not_exist(fname_default_ct_coord, '*.csv', 'Select CT coordinates table'));
    element_info    = readtable(uigetfile_if_not_exist(fname_default_element_info, '*.csv', 'Select element_info table'));
    ROIC_filename   = uigetfile_if_not_exist(fname_default_ROIC, '*.mat', 'Select ROI Centers file');
    fs_subj_dir     = uigetfile_if_not_exist(fname_default_fs_subj, [], 'Choose freesurfer subject directory (e.g. ~/freesurfer/subj1/', USE_SYM_LINK, 1);
    mr_path         = uigetfile_if_not_exist(locDirs.mr_pre_dicom, [], 'Optional: choose MR dicom directory (skip by closing dialog if you have .nii already)', 0, 1);
    ct_path         = uigetfile_if_not_exist(locDirs.ct_1_dicom, [], 'Optional: choose CT dicom directory (skip by closing dialog if you have .nii already)', 0, 1);
    
    if isempty(mr_path) || is_empty_dir(mr_path)
        uigetfile_if_not_exist(fname_default_mr_nii, '*.nii', 'Select MR nifti');
    end
    if isempty(ct_path) || is_empty_dir(ct_path)
        uigetfile_if_not_exist(fname_default_ct_nii, '*.nii', 'Select CT nifti');
    end
    
    try
        temp = load(ROIC_filename);
        std_ROICs = temp.vert_idx(1:N_ROIC);
    catch e
        fprintf('Error: could not load ROICs: %s\n', e.message);
    end
        
    localizer_get_status(subj, rootDir, 1);

    read_dicoms(locDirs.mr_pre_dicom, 'mr_pre', locDirs.mr_pre,'afni_bin',AFNI_BIN, 'freesurfer_bin',FS_BIN);
    read_dicoms(locDirs.ct_1_dicom, 'ct_implant', locDirs.ct_1);
    create_surf(subj,fullfile(locDirs.mr_pre,'mr_pre.nii'), locDirs.fs, 'freesurfer_bin',FS_BIN);
    suma(subj, locDirs.fs, 'afni_bin',AFNI_BIN, 'freesurfer_bin',FS_BIN);
    
    dural_filename = fullfile(locDirs.fs_subj,'SUMA',[hemisphere 'h.pial-outer-smoothed.gii']);
    pial_filename = fullfile(locDirs.fs_subj,'SUMA',['std.141.' hemisphere 'h.pial.gii']);
    assert(exist(dural_filename, 'file') > 0, 'File not found: %s\n', dural_filename);
    assert(exist(pial_filename, 'file') > 0, 'File not found: %s\n', pial_filename);
    dural_surf = gifti(dural_filename);
    pial_surf = gifti(pial_filename);
    
    %% anchors
    if ~exist(fullfile(locDirs.anchors, 'anchors.csv'), 'file')
        anchor_table = anchor(fullfile(locDirs.fs, subj), fullfile(locDirs.mr_pre,'mr_pre.nii'), locDirs.anchors,'slicer_binary',SLICER_BINARY,'postprocess_only',0);
        writetable(anchor_table, fullfile(locDirs.anchors, 'anchors.csv')); 
    else
        anchor_table = readtable(fullfile(locDirs.anchors, 'anchors.csv'));
    end
    
    ct_file = register_CT(fullfile(locDirs.mr_pre,'mr_pre.nii'), fullfile(locDirs.ct_1, 'ct_implant.nii'), locDirs.ct_1_xfm);
    leads = xfm_leads(ct_coords, fullfile(locDirs.ct_1_xfm, ct_file), fullfile(locDirs.ct_1, 'ct_implant+orig'), 'Voxel', locDirs.ct_1_xfm);
    proj_coords = projection(leads, dural_surf, anchor_table,element_info);
    writetable(proj_coords, fullfile(locDirs.mr_pre, 'coords.csv'));
    roi_mesh_d_lut = grow_ROIs(pial_surf, std_ROICs, R_ROI);
    lead_roi_lut = lead_to_ROI(roi_mesh_d_lut, proj_coords, pial_surf, dural_surf, R_ROI, locDirs.roi, 'contact_radius',R_CONTACT);
    writetable(lead_roi_lut, fullfile(locDirs.roi, 'lead_roi_lut.csv'));
    writetable(roi_mesh_d_lut, fullfile(locDirs.roi, 'roi_mesh_d_lut.csv'));
    
    %% Display results of projection
    f = figure;
    g = gifti(fullfile(locDirs.fs_subj,'SUMA',sprintf('std.141.%sh.pial.gii',hemisphere)));
    p = patch('Vertices',surf.vertices,'Faces',surf.faces);
    set(p,'FaceColor',0.8*[1 1 1]);
    set(p,'EdgeAlpha',0);
    camlight; axis equal; hold on
    
    xyz = proj_coords{:,{'x','y','z'}};
    hw = unique( cellfun(@(s) s(isletter(s)), proj_coords.chanName, 'uniformOutput',0), 'stable' );
    colors = hsv(numel(hw));
    for i = 1:numel(hw)
        ei = strncmp(proj_coords.chanName, hw{i}, length(hw{i}));
        plot3(xyz(ei,1),xyz(ei,2),xyz(ei,3),'ko','markerSize',20,'markerFaceColor',colors(i,:))
        plot(xyz(ei,1),xyz(ei,2),xyz(ei,3),'k');
    end

end
    
function filepath = uigetfile_if_not_exist(default_path, filespec, title, sym, isdir)

    if nargin < 3, title = ''; end
    if nargin < 4, sym = 0; end
    if nargin < 5, isdir = 0; end
    
    if exist(default_path, 'file') && ~(isdir && is_empty_dir(default_path))
        filepath = default_path;
    else
        fprintf('%s\n', title);
        if isdir
            filepath = uigetdir(default_path, title);
        else
            [filename,pathname] = uigetfile(filespec, title, fileparts(default_path));
            filepath = fullfile(pathname, filename);
        end
        if ischar(filepath) && exist(filepath, 'file') && ~strcmpi(filepath, default_path)
            if sym
                if exist(default_path, 'dir')
                    rmdir(default_path,'s');
                end
                fprintf('Linking %s <-- %s\n', filepath, default_path);
                unix(sprintf('ln -s %s %s', filepath, default_path));
            else
                fprintf('Copying %s --> %s\n', filepath, default_path);
                copyfile(filepath, default_path);
            end
        else
            filepath = [];
        end
    end
end

function empty = is_empty_dir(filepath)
    s = dir(filepath);
    charname = char(s.name);
    count = sum(charname(:,1) ~= '.');
    empty = count == 0;
end
