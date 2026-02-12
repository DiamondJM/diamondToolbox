function anchor_table = anchor(subj_dir, mr_filename, varargin)
% anchor 
%
% Usage: ANCHOR(subj_dir, ...)
%
% Inputs:
%   subj_dir        - Subject directory which must *contain* SUMA/ folder with resampled surfaces
%   mr_filename     - Path to MR .nii volume file
%   
% Optional Input:
%   working_dir     - Output/intermediate file directory. Default is current directory.
%   
%   
%
% Input key-values:
%   log_file        - path to log file. 
%   preprocess_only - If true, only prepare the files needed for slicer
%   postprocess_only- If true, only convert the fiducial markers from fcsv to csv
%   slicer_binary   - Full path to Slicer GUI executable (eg. Slicer/Contents/MacOS/Slicer
%   subject         - subject ID. By default, this value is taken to be equal to the subj_dir folder name
%   
% Output: 
%   anchor_table    - Table with chanName and x,y, and z coordinate columns for anchors. Columns are:
%       chanName      - Channel name e.g. 'G2'
%       x, y, and z   - coordinates
%       tagName       - Tagname (hardware) e.g. 'G'
%       absElmtIdx    - Index of channel within hardware e.g. 2
%       orientation   - LPI
%       subjId        - subject name
%       useAnchor     - Default 1. Can toggle to 0 to not use in projection algorithm
%       
%
% File Output:
%   pial and pial-outer_smoothed.stl files for both hemispheres
%   anchors.csv
%   anchors.mrml
%
% Description:
%   1) Prepares 3D slicer scene file for user to manually identify anchor coordinates.
%   2) opens 3D Slicer.
%   3) Converts Slicer's fcsv file to csv table of xyz-coordinates.
%   
%

% Revision History:
%   03/17 - MST

%  Copyright (C) 2017 Mike Trotta

    disp anchor
    
    % Input
    ip = inputParser;
    ip.addOptional('working_dir', pwd, @(x) ischar(x) && exist(x,'file')>0);
    ip.addParameter('preprocess_only', false);
    ip.addParameter('postprocess_only', false);
    ip.addParameter('log_file', []);
    ip.addParameter('slicer_binary', []);
    ip.addParameter('subject', []);
    ip.parse(varargin{:});
    working_dir     = ip.Results.working_dir;
    log_file        = ip.Results.log_file;
    preprocess_only = ip.Results.preprocess_only;
    postprocess_only = ip.Results.postprocess_only;
    slicer_binary   = ip.Results.slicer_binary;
    subject         = ip.Results.subject;
    
    anchor_table = [];
    files.fname_lpialgii = 'lh.pial.gii';
    files.fname_rpialgii = 'rh.pial.gii';
    files.fname_lpialstl = 'lh_pial.stl';
    files.fname_rpialstl = 'rh_pial.stl';
    files.fname_lenvgii = 'lh.pial-outer-smoothed.gii';
    files.fname_renvgii = 'rh.pial-outer-smoothed.gii';
    files.fname_lenvstl = 'lh_pial-outer-smoothed.stl';
    files.fname_renvstl = 'rh_pial-outer-smoothed.stl';
    %files.mr_pre        =  mr_filename;
    mrml_name           = 'anchors.mrml';
    scene_file = fullfile(working_dir, mrml_name);
    fiducial_name       = 'anchors.fcsv';
    suma_dir = fullfile(subj_dir, 'SUMA');
    
    %assert(exist(mr_filename,'file') > 0, 'File not found: %s\n', mr_filename);
    %copyfile(mr_filename, working_dir);
    
    % Open log
    try fd = fopen(log_file, 'a+');
    catch, fd = 1;
    end
    
    try
        % Prepare .mrml scene file (ie preprocess)
        if ~postprocess_only 
            fprintf(fd, 'Creating Slicer .mrml scene file...\n');
            
            if ~exist(suma_dir, 'dir')
                error('Required SUMA directory not found: %s\n', suma_dir);
            end
            
            % Check existance of gifti files, prepend full path
            fields = fieldnames(files);
            for i = 1 : length(fields)
                filename = files.(fields{i});
                [~,~,ext] = fileparts(filename);
                if strcmp(ext, '.gii') 
                    filename = fullfile(suma_dir, filename);
                    files.(fields{i}) = filename;
                    if ~exist(filename, 'file')
                        if ~inputYN(sprintf('Required file not found: %s\nContinue anyway? ', filename))
                            error('Required file not found: %s\n', filename);
                        end
                    end
                    
                elseif strcmp(ext, '.stl')
                    filename = fullfile(working_dir, filename);
                    files.(fields{i}) = filename;

                end
                
            end
            
            exist_all_stl = count_filetype(working_dir, 'stl') == 4;
            
            % Run File conversion command
            if ~exist_all_stl
                command = sprintf(['ConvertSurface -i_gii %s -o_stl %s -orient_out "LPI"; ' ...
                        'ConvertSurface -i_gii %s -o_stl %s -orient_out "LPI"; ' ...
                        'ConvertSurface -i_gii %s -o_stl %s -orient_out "LPI"; ' ...
                        'ConvertSurface -i_gii %s -o_stl %s -orient_out "LPI"; '], ...
                        files.fname_lpialgii, files.fname_lpialstl, ...
                        files.fname_rpialgii, files.fname_rpialstl, ...
                        files.fname_lenvgii, files.fname_lenvstl, ...
                        files.fname_renvgii, files.fname_renvstl);

                if ~isempty(log_file), status = bashlog2(command, log_file);
                else, [status,txt] = unix(command);
                end

                if status ~= 0
                    fprintf('%s\n',txt);
                    error('Creating stl files from gifti failed'); 
                end
            end
            
            % Create scene files from template
            if ~count_filetype(working_dir, 'mrml')
                mrml_path = which('template_anchors.mrml');
                copyfile(mrml_path, scene_file); 
            end
            
            if ~count_filetype(working_dir, 'fcsv')
                fcsv_path = which('template_anchors.fcsv');
                copyfile(fcsv_path, fullfile(working_dir, fiducial_name));
            end    
            
        end
            
        

        % Launch slicer for user
        if ~postprocess_only || ~preprocess_only
            [result,~] = unix('type Slicer');
            is_slicer_on_path = (0 == result);
            if is_slicer_on_path
                script = 'Slicer'; 
            elseif ~isempty(slicer_binary) && unix(['type ' slicer_binary]) == 0
                script = slicer_binary;
            else
                script = [];
            end
            if ~isempty(script)
                
                fprintf('\n\n-------------------------------------------------------------------\n');
                fprintf('Launching external Slicer application. Please wait..\n');
                fprintf('Matlab script is blocked and will continue upon closing application.\n');
                
                command = sprintf('%s -i %s &> /dev/null', script, scene_file);
                
                fprintf('\nInstructions:\n');
                fprintf('\tPlace fiducial markers now. Marker names must correspond to channel names.\n');
                fprintf('\tSave fiducial file as "anchors.fcsv" in the working_dir (%s)\n', working_dir);
                fprintf('\tThen *close Slicer* and continue with this script.\n');
                fprintf('-------------------------------------------------------------------\n');
            
                unix(command);
                
                if ~inputYN('Did you finished marking all anchors and save an anchor.fcsv file?')
                    fprintf('You may re-call this function with postprocess_only=1 when you have finished marking\n');
                    return;
                end
                
            else
                fprintf('Could not find Slicer on the path. Launch Matlab from terminal to inherit PATH, or\n');
                fprintf('pass in slicer_binary argument. Alternatively, launch application yourself now...\n');
                pause(5);
            end
                
            % Instruction to user
            
        end

        % Convert fcsv --> table (ie postprocess)
        if ~preprocess_only
            fprintf(fd, 'Converting fcsv fiducial coordinates to csv table\n');
            fcsv = fullfile(working_dir, 'anchors.fcsv');
            if ~exist(fcsv, 'file')
                error('Needed file not found: %s\n', fcsv);
            end
            
            % Infer from subj_dir
            if isempty(subject)
                [~,subject] = fileparts(subj_dir);
            end
            
            anchor_table = slicer_read_fiducials(fcsv);
            if isempty(anchor_table)
                warning('No anchors found in table (%s)', fcsv)
            end
            
            anchor_table.subjId = repmat({subject},height(anchor_table),1);
            anchor_table.useAnchor = ones(height(anchor_table),1);
            
        
        end
        
    catch e
        if fd > 1, fclose(fd); end
        rethrow(e);
    end
end
