function mr_coords = xfm_leads(ct_coords, xfm_file, ct_file, coord_format, varargin)
% xfm_leads 
%
% Usage: mr_coords = xfm_leads(ct_coords, xfm_file, ct_file, ...)
%
% Inputs:
%   ct_coords      - table with chanName, x, y, and z columns representing original registered CT electrode coordinates
%   xfm_file       - transformation matrix output from CT to MR registration
%   ct_file        - path to original AFNI post-implant CT file (e.g. /path/to/ct_implant+orig)
%   coord_fmt      - format of ct_coords. Either 'voxel' or 'RAI' (voxels means coordinates in RAI-ijk format, RAI means millimeter xyz format)
%
% Optional Input:
%   working_dir     - directory where intermediate files and output files will be saved
%   
%
% Input key-values:
%   log_file        - path to log file. 
%   
% Output: 
%   mr_coords       - A table like ct_coords but with the x,y,z columns transformed to MR space via xfm_file
%
% File Output:
%   vox2mm.1D       - temporary file encoding the CT's voxel (ijk) to real millimeter xyz transformation
%   vox.txt         - temporary file voxel coordinates of original CT
%   mm.txt          - temporary file millimeter coordinates of original CT
%   elec.gii        - temporary file giving electrodes as a surface, transformed to MR via given transformation
%   
% Description:
%   Applies the transform representing CT to MR space to CT-registered coordinates
%   
%
% Revision History:
%   03/17 - MST

%  Copyright (C) 2017 Mike Trotta

    fname_surf = 'elec.gii';
    
    % -------------
    %%--- Setup ---
    % -------------
    
    % Input
    ip = inputParser;
    ip.addOptional('working_dir', pwd, @ischar);
    ip.addParameter('log_file', []);
    ip.parse(varargin{:});
    working_dir = ip.Results.working_dir;
    log_file    = ip.Results.log_file;
    
    % Check input
    required_cols = {'chanName','x','y','z'};
    assert(istable(ct_coords) && length(required_cols) == length(intersect(required_cols, ct_coords.Properties.VariableNames)),...
        'ct_coords must be a table with chanName, x, y, and z columns');

    
    assert(exist(xfm_file, 'file') > 0, 'xfm_file not found: %s\n', xfm_file);
    assert(exist(sprintf('%s.BRIK', ct_file), 'file') > 0, 'CT file not found: %s.BRIK\n', ct_file);
    
    switch upper(coord_format)
        case 'VOXEL',   fname_vox = 'vox.txt';  is_vox = 1;
        case 'RAI',     fname_vox = 'mm.txt';   is_vox = 0;
        otherwise, error('Unkown coord_format %s. must be Voxel or RAI', coord_format);
    end    
    
    % ---------------------------------------
    % Apply transformation to CT coordinates
    % ---------------------------------------
    
    % TODO: here is where we could switch the implementation to use vecWarp directly instead of this gifti surface
    
    % Write the voxel coordinates in a file
    xyz = ct_coords(:,{'x','y','z'});
    writetable(xyz, fullfile(working_dir, fname_vox),...
        'WriteVariableNames',false,...
        'Delimiter','space');
    
    % Call the shell script
    dname_loc = fileparts(fileparts(mfilename('fullpath'))); % up two directories
    script = fullfile(dname_loc, 'shell_scripts', 'xfm_leads.sh');
    
    if is_vox
        command = sprintf('%s %s %s %s %s %s', script, working_dir, ct_file, xfm_file, fname_vox);
    else
        command = sprintf('%s %s %s %s %s', script, working_dir, ct_file, xfm_file);
    end

    if ~isempty(log_file), status = bashlog(command, log_file);
    else, status = unix(command);
    end

    if status ~= 0, error('Transformation script failed'); end
    
    % Extract coordinates from gifti, store in table
    mr_xyz = nimg_gifti2coords(fullfile(working_dir,fname_surf));
    mr_coords = ct_coords;
    mr_coords{:, {'x','y','z'}} = mr_xyz{:,:};
    
end
