function read_dicoms(d_dicom, name, varargin)
% read_dicoms converts a directory of dicoms to NIFTI and AFNI formats
%
% Usage: READ_DICOMS(d_dicom, name, working_dir, ...)
%
% Inputs:
%   d_dicom         - directory containing dicoms
%   name            - name of output image (e.g. 'subj' results in 'subj.nii' output)
%   
% Optional Input:
%   working_dir     - Where files should be put. Default: current directory.
%   
%
% Input key-values:
%   log_file        - path to log file. 
%   freesurfer_bin  - path to freesurfer bin. Needed if matlab not started from terminal.
%   afni_bin        - path to freesurfer bin. Needed if matlab not started from terminal.  
%
% Output:   
%   None
%
% File Output:
%   name.nii
%   name.BRIK
%   name.HEAD
%
% Description:
%   Uses Freesurfer's mri_convert and AFNI's 3dCopy to convert dicoms to nifti and AFNI formats
%   
%

% Revision History:
%   03/17 - MST

%  Copyright (C) 2017 Mike Trotta
    
    ip = inputParser;
    ip.addOptional('working_dir', pwd, @ischar);
    ip.addParameter('log_file', []);
    ip.addParameter('freesurfer_bin', []);
    ip.addParameter('afni_bin', []);
    ip.parse(varargin{:});
    working_dir = ip.Results.working_dir;
    log_file = ip.Results.log_file;
    afni_bin = ip.Results.afni_bin;
    fs_bin = ip.Results.freesurfer_bin;
    
    % function counts the files with given extension in a directory
    
    if ~exist(working_dir, 'dir'), error('Not found: %s', working_dir); end
    if ~exist(d_dicom,'dir'), error('Not found: %s', d_dicom); end
    
    nifti_exists = count_filetype(working_dir, 'nii') > 0;
    dicoms_exist = count_filetype(d_dicom, 'dcm') > 0;
    afni_exists  = count_filetype(working_dir, 'BRIK') > 0;
    
    if ~dicoms_exist && ~nifti_exists
        error('read_dicoms: No dicoms (.dcm) or nifti (.nii) found');
    end
    
    dname_this = fileparts(fileparts(mfilename('fullpath'))); % up 2 dirs
    fname_script = fullfile(dname_this, 'shell_scripts', 'read_dicoms.tcsh');
    assert(exist(fname_script,'file') > 0, 'File not found: %s', fname_script);

    % Check path for binaries
    path_os = getenv('PATH');
    if ~isempty(fs_bin), path_os = strcat(':',path_os,':',fs_bin,':'); end
    if ~isempty(afni_bin), path_os = strcat(':',path_os,':',afni_bin,':'); end
    if ~isempty(fs_bin)
        fs_text = 'freesurfer';
        ndx = strfind(fs_bin, fs_text) + length(fs_text);
        fs_home = fileparts(fs_bin(1:ndx));
        setenv('FREESURFER_HOME', fs_home);
    end
    setenv('PATH',path_os);
    
    if unix('type afni') ~= 0 || unix('type freesurfer') ~= 0
        error(['ERROR: Cant find the afni or freesurfer command from MATLAB'...
            '. Start MATLAB from terminal or pass in bin path']);
    end
    
    % Execute shell script
    if nifti_exists && ~dicoms_exist && ~afni_exists
        % dicoms are not actually necessary. do .nii -> brik/head directly
        command = sprintf('3dcopy %s.nii %s', fullfile(working_dir, name), fullfile(working_dir, name));
    else
        command = sprintf('%s %s "%s" "%s"', fname_script, name, d_dicom, working_dir);
    end
    fprintf('%s\n',command);
    
    if ~isempty(log_file)
        status = bashlog2(command, log_file);
    else
        [status,msg] = unix(command);
        fprintf('%s',msg);
    end
    
    if status ~= 0
        error('read_dicoms failed: %s',msg);
    end
    
        
end

