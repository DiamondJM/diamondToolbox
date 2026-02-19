function status = localizer_get_status(subj, rootEEGdir, report)
% status = localizer_get_status(subj, rootEEGdir, report)
%
%   Inputs
%       subj       - code name of subject
%       rootEEGdir - rootEEGdir should be folder containing subject folders
%       report     - Optional, default 0. If 1, print a basic report of number of files
%
% This file is intended to be a template for your own localization pipeline. It serves to
% return (and if report=1) display the number of files in the subject's directory structure
% matching certain file extensions. It is useful to make sure the required input files are 
% in place.
%
%

%  Copyright (C) 2017 Mike Trotta
% 

    if nargin < 3, report = 0; end
    
    locDirs = localizer_create_directories(subj, rootEEGdir);
    
    % check for MR and CT dicom/niftis
    status.n_mr_dcm = count_filetype(locDirs.mr_pre_dicom, {'dcm','dicom'});
    status.n_mr_dcm2 = count_filetype(locDirs.mr_t2_dicom, {'dcm','dicom'});
    status.n_mr_nii = count_filetype(locDirs.mr_pre, {'nii'});
    status.n_mr_nii2 = count_filetype(locDirs.mr_t2, {'nii'});
    status.n_mr_afni= count_filetype(locDirs.mr_pre, {'brik'});
    status.n_ct_dcm = count_filetype(locDirs.ct_1_dicom, {'dcm','dicom'});
    status.n_ct_nii = count_filetype(locDirs.ct_1, {'nii'});
    status.n_ct_afni= count_filetype(locDirs.ct_1, {'brik'});
    if isfield(locDirs, 'ct_t')
        status.n_ct_dcm2 = count_filetype(locDirs.ct_2_dicom, {'dcm','dicom'});
        status.n_ct_nii2 = count_filetype(locDirs.ct_2, {'nii'});
        status.n_ct_afni2= count_filetype(locDirs.ct_2, {'brik'});
    end
    
    % check freesurfer surfaces and outer smoothed surfaces + suma
    d_surf   = fullfile(locDirs.fs_subj, 'surf');
    status.n_pial   = count_filetype(d_surf, {'pial'});
    status.n_pial_smooth = count_filetype(d_surf, {'pial-outer-smoothed'});
    
    % anchors
    status.n_stl    = count_filetype(locDirs.anchors, {'stl'});
    
    if report
        fprintf('---------------------------------\n');
        fprintf('| Localizer files - %6s      |\n', subj);
        report_line('MR T1 dicoms', status.n_mr_dcm);
        report_line('MR T1 nifti', status.n_mr_nii);
        report_line('MR T2 dicoms', status.n_mr_dcm2);
        report_line('MR T2 nifti', status.n_mr_nii2);
        report_line('MR afni', status.n_mr_afni);
        report_line('CT 1 dicom', status.n_ct_dcm);
        report_line('CT 1 nifti', status.n_ct_nii);
        report_line('CT 1 afni', status.n_ct_afni);
        if isfield(locDirs, 'ct_t')
            report_line('CT 2 dicom', status.n_ct_dcm2);
            report_line('CT 2 nifti', status.n_ct_nii2);
            report_line('CT 2 afni', status.n_ct_afni2);
        end
        report_line('pial surf', status.n_pial);
        report_line('pial-smooth surf', status.n_pial_smooth);
        report_line('anchor stl', status.n_stl);
        fprintf('---------------------------------\n');
    end
end



function report_line(filetype, count)
    if count == 0
        fprintf('| %-18s: [MISSING] |\n',filetype)
    else
        fprintf('| %-18s: [  %2d   ] |\n',filetype, count)
    end
end