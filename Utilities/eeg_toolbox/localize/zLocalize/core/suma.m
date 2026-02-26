function suma(subject, varargin)
% suma 
%
% Usage: suma(subject, working_dir, ...)
%
% Inputs:
%   subject         - name of subject
%   
% Optional Input:
%   working_dir     - Where files should be put. Default: current directory. This should
%                     be the path to a freesurfer directory with the subject name or the directory 
%                     which contains this subject's freesurfer folder.
%   
%
% Input key-values:
%   log_file        - path to log file. 
%   afni_bin        - path to afni bin. Needed if matlab not started from terminal.
%   freesurfer_bin  - path to freesurfer bin. Needed if matlab not started from terminal.
%   rerun           - if 1, rerun
%
% Output: 
%   None
%
% File Output:
%   Suma files in subject's freesurfer folder
%
% Description:
%   Runs SUMA's @SUMA_Make_Spec_FS which resamples surface mesh to standard
%

% Revision History:
%   03/17 - MST

%  Copyright (C) 2017 Mike Trotta

    % Input
    ip = inputParser;
    ip.addOptional('working_dir', pwd, @ischar);
    ip.addParameter('log_file', []);
    ip.addParameter('afni_bin', []);
    ip.addParameter('freesurfer_bin', []);
    ip.addParameter('rerun', []);

    ip.parse(varargin{:});
    working_dir = ip.Results.working_dir;
    log_file    = ip.Results.log_file;
    afni_bin    = ip.Results.afni_bin;
    fs_bin      = ip.Results.freesurfer_bin;
    rerun      = ip.Results.rerun;
    
    
    % Make working_dir
    if ~exist(working_dir, 'dir'), mkdir(working_dir); end
    
    [d_parent, d] = fileparts(working_dir); % d is just name of folder
    if strcmp(d, subject)
        % working_dir is the subject folder
        subj_root = d_parent;
    else
        % working_dir will contain the subject folder
        subj_root = working_dir;
    end
    
    subj_folder = fullfile(subj_root, subject);
    if ~exist(subj_root, 'dir'), mkdir(subj_root); end
    if ~exist(subj_folder, 'dir'), mkdir(subj_folder); end
    
    % Check path for binaries
    if ~check_for_system_bin('suma')
        addpath_system(afni_bin);
    end
    if ~check_for_system_bin('mri_convert')
        addpath_system(fs_bin);
    end

    if ~check_for_system_bin('suma')
        error(['ERROR: Cant find the suma command from MATLAB'...
            '. Start MATLAB from terminal or pass in afni_bin']);
    end
    
    % clear files necessary to rerun
    if rerun
        delete(fullfile(subj_folder, 'SUMA', sprintf('%s_lh.spec',subject)));
        delete(fullfile(subj_folder, 'SUMA', sprintf('%s_rh.spec',subject)));
        delete(fullfile(subj_folder, 'SUMA', sprintf('%s_SurfVol.nii',subject)));
    end
    
    % Run 
    command = sprintf('SUBJECTS_DIR=%s; @SUMA_Make_Spec_FS -sid %s -fspath %s -GIFTI -ld 141', ...
        subj_root, subject, subj_folder);
    
        
    if ~isempty(log_file), status = bashlog2(command, log_file);
    else, status = unix(command);
    end

    if status ~= 0, error('SUMA command failed'); end
    
    
end

function addpath_system(p)
    setenv('PATH',strjoin({getenv('PATH'), p}, ':'));
end

function is_on_path = check_for_system_bin(s)
    is_on_path = (0 == system(sprintf('which %s', s)));
end