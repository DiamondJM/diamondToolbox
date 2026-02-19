function create_surf(subject, mr_nifti, varargin)
% create_surf creates pial and envelope surfaces from an MRI volume
%
% Usage: CREATE_SURF(subject, mr_nifti, working_dir, ...)
%
% Inputs:
%   subject         - name of subject
%   mr_nifti        - filepath to MR nitfti to use for surface. MR should be T1 MPRAGE.
%   
% Optional Input:
%   working_dir     - Where files should be put. Default: current directory. This should
%                     be the path to a freesurfer directory with the subject name or the directory 
%                     which contains this subject's freesurfer folder.
%   
%
% Input key-values:
%   log_file        - path to log file. 
%   freesurfer_bin  - path to freesurfer bin. Needed if matlab not started from terminal.
%   envelope_only   - If true, skip making the surface with recon-all (useful if surface already created)
%   T2_nifti        - Pass in the path to a T2 nifti, which improves surface reconstruction
%   
% Output: 
%
% File Output:
%   Directory with subject name and containing freesurfer files
%   Cortical envelope output as pial-outer-smoothed surface
%   
%
% Description:
%   Creation of pial:
%       Freesurfer's recon-all
%
%   Creation of pial envelope:
%       Execute freesurfer local-gi procedure. Calls subfunctions explicitly,
%       thereby saving computaitonal time.
%       See https://surfer.nmr.mgh.harvard.edu/fswiki/LGI for a description.
%       Commands are carried out in the SUBJ/surf freesurfer directory.
%
%  Copyright (C) 2017 Mike Trotta

% Revision History:
%   03/17 - MST
%   06/10/2020 SJ: added the following to allow for freesurfer update 7.1.0 or future updates where 
%                  freesurfer is within another file with the version number
%

    CLOSE_SPHERE_VOLUME = 12; % original default 15
    DILATION_RADIUS = 2; % note we must dilate because smoothing erodes

    % Input
    ip = inputParser;
    ip.addOptional('working_dir', pwd, @ischar);
    ip.addParameter('log_file', []);
    ip.addParameter('freesurfer_bin', []);
    ip.addParameter('envelope_only', 0);
    ip.addParameter('T2_nifti', []);
    ip.parse(varargin{:});
    working_dir = ip.Results.working_dir;
    log_file = ip.Results.log_file;
    fs_bin = ip.Results.freesurfer_bin;
    skip_recon = ip.Results.envelope_only;
    T2_nifti = ip.Results.T2_nifti;
    
    if ~exist(mr_nifti, 'file'), error('Not found: %s', mr_nifti); end
    if ~isempty(T2_nifti) && ~exist(T2_nifti, 'file'), error('Not found: %s', T2_nifti); end
    if ~exist(working_dir, 'dir'), mkdir(working_dir); end
    
    % Var / directory setup
    [d_parent,d] = fileparts(working_dir); % d is just name of folder
    if strcmp(d, subject)
        % working_dir is the subject folder
        subj_root = d_parent;
    else
        % working_dir will contain the subject folder
        subj_root = working_dir;
    end
    
    dname_this = fileparts(mfilename('fullpath'));
    dname_loc  = fileparts(dname_this);
    
    subj_folder = fullfile(subj_root, subject);
    if ~exist(subj_root, 'dir'), mkdir(subj_root); end
    if exist(subj_folder, 'dir')
        if numel(dir(subj_folder)) <= 2 % (. and .. always present)
            rmdir(subj_folder); 
        else
            fprintf('subject folder exists and is non-empty: %s\nTherefore skipping recon-all\n', subj_folder);
            skip_recon = 1;
        end
    end
    
    % Check path for binaries
    path_os = getenv('PATH');
    if ~isempty(fs_bin)
        path_os = strcat(':',path_os,':',fs_bin); 
        fs_text = 'freesurfer';
        %SJ: added the following to allow for freesurfer update 7.1.0 or future updates where freesurfer
        %is within another file with the version number
        fsvernum = char(regexp(fs_bin,['(?<=' fs_text filesep ').*(?=' filesep 'bin)'],'match'));
        if ~isempty(fsvernum)
            fs_text = [fs_text filesep fsvernum];
        end
            
        ndx = strfind(fs_bin, fs_text) + length(fs_text);
        fs_home = fileparts(fs_bin(1:ndx));
        setenv('FREESURFER_HOME', fs_home);
    end
    setenv('PATH',path_os);
    
    if  unix('type freesurfer') ~= 0
        error(['ERROR: Cant find the freesurfer command from MATLAB'...
            '. Start MATLAB from terminal or pass in bin path']);
    end
    
    % ---------------------------------------
    % Run command for recon-all, creating pial surface
    % ---------------------------------------
    if ~skip_recon
        setenv('SUBJECTS_DIR', subj_root);
        command = sprintf('recon-all -sd %s -subjid %s -all -i %s',...
            subj_root, subject, mr_nifti);
        
        if ~isempty(T2_nifti)
            % Align T2 to T1
            fname_script    = fullfile(dname_loc, 'shell_scripts', 'proc_t2_mr.sh');
            [~,temp]        = fileparts(T2_nifti);
            T2_aligned_nii  = sprintf('%s_aligned.nii', temp);
            preproc_command = sprintf('%s %s %s %s %s', fname_script, mr_nifti, T2_nifti, T2_aligned_nii, working_dir);
            command         = sprintf('%s; %s -T2 %s -T2pial', preproc_command, command, T2_aligned_nii);
        end
        fprintf('%s\n', command);
        
        if ~isempty(log_file), status = bashlog2(command, log_file);
        else, [status,txt] = unix(command);
        end

        if status ~= 0
            fprintf('%s',txt);
            error('create_surf recon-all failed'); 
        end
        
        %% Lobes
        cur_pwd = pwd;
        
        cd(fullfile(subj_root, subject));
        command = sprintf([ 'mri_annotation2label --s %s --hemi rh --sd %s --ctab lobesLUT_rh.txt --lobesStrictPHCG rh.aparc.lobesStrict; ' ...
                            'mri_annotation2label --s %s --hemi lh --sd %s --ctab lobesLUT_lh.txt --lobesStrictPHCG lh.aparc.lobesStrict; '], ...
                            subject, subj_root, ...
                            subject, subj_root);
        system(command);
        command = sprintf([ 'mri_aparc2aseg --s %s  --labelwm --annot  aparc.lobesStrict; ' ...
                            'mri_convert mri/aparc.lobesStrict+aseg.mgh mri/aparc.lobesStrict+aseg.nii'] , ...
            subject);
        system(command);
        cd(cur_pwd);
                
        
        
    end
    
    % ---------------------------------------
    % Create cortical envelope using local_gi
    % ---------------------------------------
    
    % Check matlab path
    if isempty(which('make_outer_surface')) && ~isempty(fs_bin)
        fs_text = 'freesurfer';
        %SJ: added the following to allow for freesurfer update 7.1.0 or future updates where freesurfer
        %is within another file with the version number
        fsvernum = char(regexp(fs_bin,['(?<=' fs_text filesep ').*(?=' filesep 'bin)'],'match'));
        if ~isempty(fsvernum)
            fs_text = [fs_text filesep fsvernum];
        end
        ndx = strfind(fs_bin, fs_text) + length(fs_text);
        fs_home = fileparts(fs_bin(1:ndx));
        setenv('FREESURFER_HOME', fs_home);
        addpath(fullfile(fs_home, 'matlab'));
    end
    
    cur_pwd = pwd;
    
    % Two halves of the script
    fname_script_1 = fullfile(dname_loc, 'shell_scripts', 'make_envelope_1.tcsh');
    fname_script_2 = fullfile(dname_loc, 'shell_scripts', 'make_envelope_2.tcsh');
    assert(exist(fname_script_1,'file') > 0, 'File not found: %s', fname_script_1);
    assert(exist(fname_script_2,'file') > 0, 'File not found: %s', fname_script_2);
    
    % Execute script for each hemisphere separately
    surfs = {'lh.pial','rh.pial'};
    for i = 1 : 2, surf = surfs{i};
        cd(fullfile(subj_folder, 'surf'));
        command = sprintf('%s --i %s', fname_script_1, surf);
        fprintf('%s\n', command);
        
        % script part 1
        if ~isempty(log_file), status = bashlog2(command, log_file);
        else, [status,txt] = unix(command);
        end

        temp_dir = fullfile(pwd, sprintf('tmp-mris_compute_lgi-%s', surf)); % script creates this
        volume = fullfile(temp_dir, sprintf('%s.filled.mgz', surf));
        envelope = fullfile(temp_dir, sprintf('%s-outer', surf));
        if status ~= 0 || ~exist(temp_dir,'dir')
            cd(cur_pwd);
            fprintf('%s',txt);
            error('Create envelope script failed'); 
        end
        
        % matlab part
        
        try make_outer_surface2(volume, CLOSE_SPHERE_VOLUME, envelope, DILATION_RADIUS);
        catch e
            cd(cur_pwd);
            rethrow(e);
        end
        
        % script part 2
        command = sprintf('%s --i %s', fname_script_2, surf);
        fprintf('%s\n', command);
        if ~isempty(log_file), status = bashlog2(command, log_file);
        else, [status,txt] = unix(command);
        end
        
        if status ~= 0
            cd(cur_pwd);
            fprintf('%s',txt);
            error('Create envelope script failed'); 
        end
    end
    
    cd(cur_pwd);
end
