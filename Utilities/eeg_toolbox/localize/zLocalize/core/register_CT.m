function [fname_xfm, fname_ct_xfm] = register_CT(mr, ct, working_dir, varargin)
% register_CT aligns CT to MR and saves resulting transform
%
% USAGE: [fname_xfm, fname_ct_xfm] = register_CT(mr, ct, working_dir, ...)
%
% INPUT:
%   mr          - path to MRI nifti. Should be pre-op, T1 MPRAGE.
%   ct          - path to CT nifti (post-op).
%   
% OPTIONAL INPUT:
%   working_dir - path to directory to put output and temporary files. Default: current directory.
%
% INPUT Key-Value:
%   ct_reorient - Reorient CT to given orientation.
%                 Default: 'RAI'
%   xfm_cost    - Cost function for align_epi_anat.py. 
%                 Default: 'nmi' (mutual information)
%   remove_air  - If true, run @NoisySkullStrip on CT. Useful for noisy CT's. See Afni help.
%                 Default: false.
%   invert      - If true, invert MRI intensity (pair this with lpa cost)
%                 Default: false.
%   log_file    - path to log file. 
%   xfm_override - If passed, skip finding a transform and use this xfm file instead
%
% OUTPUT:
%   fname_xfm    - filename of CT-->MR transformation
%   fname_ct_xfm - filename of CT with transform applied (this is an AFNI brik/head fileaname)
%
% DETAILS:
%   The CT is aligned to a deobliqued version of the MR.
%   
%   MRI Processing (via AFNI functions): 
%       - deoblique MR
%       - 
%   CT Processing (via AFNI functions):
%       - reorient to given orientation (useful to match orientation of electrode coordinates)
%       - align center to MR center (rigid transform)
%       - resample to match MR
%       - remove air voxels (@NoisySkullStrip)
%   
% See also: proc_t1_mr.sh, proc_ct.sh
%

% TROUBLESHOOTING
%   Error: Cannot open dataset ??
%   Check the logged output; afni often tries to remove and rename a dataset to itself
%   look for: ** ERROR: Could not rename ....
%   cd to the working dir and rerun the most recent 3dAllineate -base... command form the log

% REVISION HISTORY:
%   03/17 MST - Created

%  Copyright (C) 2017 Mike Trotta

    XFM_NAME = 'transform.1D';
    
    % -------------
    %%--- Setup ---
    % -------------
    
    % Parameters
    ip = inputParser;
    ip.addParameter('ct_reorient','RAI');
    ip.addParameter('xfm_cost', 'nmi');
    ip.addParameter('remove_air', false);
    ip.addParameter('invert', false);
    ip.addParameter('log_file', []);
    ip.addParameter('xfm_override', []);
    ip.parse(varargin{:});
    ct_reorient         = ip.Results.ct_reorient;
    xfm_cost            = ip.Results.xfm_cost;
    remove_air          = ip.Results.remove_air;
    invert              = ip.Results.invert;
    log_file            = ip.Results.log_file;
    fname_xfm_override  = ip.Results.xfm_override;
    
    assert(nargin >= 2);
    if nargin == 2
        working_dir = pwd;
    end
    
    % Variables
    dname_loc = fileparts(fileparts(mfilename('fullpath'))); % up 2 directories
    fname_mr_proc = fullfile(dname_loc, 'shell_scripts', 'proc_t1_mr.sh');
    fname_ct_proc = fullfile(dname_loc, 'shell_scripts', 'proc_ct.sh');
    skip_xfm = false;
    
    % Check MR nifti file
    if exist(mr, 'file') && strfind(mr, '.nii') > 0
        [~,mr_name] = fileparts(mr);
        if isempty(dir(fullfile(working_dir, sprintf('%s.nii',mr_name))))
            copyfile(mr, working_dir);
        end
    else
        error('MRI nifti file not found: %s', mr);
    end
    
    % Check CT nifti file
    if exist(ct, 'file') && strfind(ct, '.nii') > 0
        [~,ct_name] = fileparts(ct);
        if isempty(dir(fullfile(working_dir, sprintf('%s.nii',ct_name))))
            copyfile(ct, working_dir);
        end
    else
        error('CT nifti file not found: %s', ct);
    end
    
    % Check xfm file if given
    if ~isempty(fname_xfm_override)
        if ~exist(fname_xfm_override, 'file')
            error('Given xfm file does not exist: %s', fname_xfm_override);
        end
        skip_xfm = true;
        fname_xfm = fname_xfm_override;
    
    end
    
    % -------------------------
    %%--- Transform section ---
    % -------------------------
    if ~skip_xfm 
        
        pwd_back = pwd;
        cd(working_dir);
        
        % Get assumed output filenames
        mr_name_deob = strcat(mr_name, '_do');
        if remove_air
            dset_ct = strcat(ct_name, '_sh2', mr_name_deob, '.rs.ma');
        else
            dset_ct = strcat(ct_name, '_sh2', mr_name_deob, '.rs');        
        end
        suffix = strcat('_XFMTO_', xfm_cost, '_', mr_name_deob);
        root_xfm = strcat(dset_ct, suffix);
        dset_ct = [dset_ct '+orig'];
        affine_xfm = strcat(root_xfm, '_mat.aff12.1D');
        fname_xfm = XFM_NAME;
        fname_ct_xfm = strcat(root_xfm, '+orig');
        if exist(fname_xfm, 'file')
            return;
        end
        
        
        % MR Processing
        command = sprintf('%s %s %s %d', fname_mr_proc, mr_name, working_dir, invert);
        if ~isempty(log_file)
            status = bashlog2(command, log_file, true);
        else
            [status] = unix(command);
        end

        % Get assumed output filenames
        if invert
            dset_mr = strcat(mr_name, '_do_amd_inv+orig');
        else
            dset_mr = strcat(mr_name, '_do_amd.ns+orig');
        end
        
        
        if isempty(dir([mr_name_deob '*']))
            error('MRI processing failed');
        end

        % CT Processing 
        if isempty(dir([dset_ct '*']))
            command = sprintf('%s %s %s %s %s', fname_ct_proc, ct_name, mr_name_deob, working_dir, ct_reorient);
            if ~isempty(log_file)
                status = bashlog2(command, log_file);
            else
                [status] = unix(command);
            end
        end

        
        
        if isempty(dir([dset_ct '*']))
            error('CT processing failed');
        end

        % Affine transform
        if ~exist(fullfile(working_dir, affine_xfm), 'file')

            strip_flag = 'None';

            command = ['align_epi_anat.py -dset1 ' dset_mr  ...
                       ' -dset2 ' dset_ct  ...
                       ' -dset2to1 '  ...
                       ' -dset1_strip ' strip_flag ...
                       ' -dset2_strip ' strip_flag ...
                       ' -suffix ' suffix  ... 
                       ' -output_dir ' working_dir  ...
                       ' -feature_size 1 ' ...
                       ' -cost ' xfm_cost  ...
                       ' -giant_move' ...
                       ' -rigid_body' ...
                       ' -overwrite'];
            disp(command);
            if ~isempty(log_file)
                bashlog2(command, log_file);
            else
                [status] = unix(command);
            end
        end
        
        % Combine affine, rigid transform to single transform named fname_xfm ("transform.1D")
        command = sprintf('cat_matvec -ONELINE %s %s_shft.1D > %s', affine_xfm, ct_name, fname_xfm);
        if ~isempty(log_file)
            bashlog2(command, log_file);
        else
            [status] = unix(command);
        end
        
        cd(pwd_back);

        
    end % End transform section
    
    
end
