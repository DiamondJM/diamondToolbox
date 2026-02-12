% copySubsOneTask
% 
% Grab and copy data for multiple subjects, but only the behavioral/eeg data associated with a
% specific task
% 
% Files/directories that are copied:
% - Any behavioral/your_task folder (including postOp and preOp)
% - any eeg/sessions_with_your_task (including noreref, processed, and
%   processedBP)
% - any other folders/files in the subject directory that are NOT a 'micro' folder or 'raw' folder.
% 
% This script also generates a copy report (copy_info.txt) that specifies
% what was copied (useful if you can't remember or aren't sure what was
% deleted in the future, or you don't know if everything copied while away
% from your computer)
% 
% I find it most helpful to run these as 3 separate sections:
%   1. STEP 1: FILL IN PARAMETERS
%       This section, which requires user input to specify parameters
%   2. STEP 2: Generate list of files/directories to copy
%       Section to generate a list of files/directories to copy. It is good
%       to check the list before giving the go-ahead for section 3 to start
%       copying
%   3. STEP 3: COPY
%       Section to actually copy the data specified in the list from the
%       previous section
% 
% Created by Samantha Jackson
% 
% 9/2019 Created by SJ to copy data for someone
% 10/2019 Modified by SJ to allow multiple sessions per task, clean to commit to eeg_toolbox
% 

% PLEASE FILL OUT THE FOLLOWING BEFORE RUNNING:

%% STEP 1: FILL IN PARAMETERS

sublist = {'NIH070'}; % Subjects you want to copy
%sublist = {'NIH043', 'NIH048', 'NIH049', 'NIH050', 'NIH051', 'NIH052', 'NIH058'}; 
subdir = '/Volumes/EEG_56STAGE-1/stage/working'; % Subject directory you want to copy from
%subdir = '/Volumes/Shares/FRNU/data/eeg';

outputdir = '/Volumes/Seagate Backup Plus Drive/forShawn_70'; % Directory you want to copy to
foldnames2 = {'behavioral', 'behavioral_postOp', 'behavioral_preOp', ...
                'eeg.noreref', 'eeg.processed', 'eeg.processedBP'}; % Folders you want to copy from
Task2Copy = 'namingTask'; % Name of the task you want copied

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2: Generate list of files/directories to copy
barray = cell(numel(sublist),3);

for subn = 1:numel(sublist) % Iterate through subjects
    sub = sublist{subn};
    fprintf('%s\n',[sub ': ']);
    % Read alignment summary to find which sessions your task is part of
    alignsum = [subdir filesep sub filesep 'behavioral' filesep 'alignmentSummary.txt'];
    alignarray = fileread(alignsum); 
    % match the expression:
    %     Find where it begings with '>>' and start matching once you hit
    %     the session number (\d{6}_\d{4}), then look ahead for
    %     'namingTask/session_'
    sess_list = regexp(alignarray,['(?<=>>[^>>]*\s)\d{6}_\d{4}(?=[^>>]*' Task2Copy '\/session_)'],'match');
    
    % Automatically first add other non-behavioral/eeg files/directories
    % like docs and tal folder:
    if isempty(sess_list)
        fprintf('%s\n',['No occurrences of ' Task2Copy 'for Subject ' sub '. Skipping to next Subject.']);
        continue;    
    else
        cd([subdir filesep sub])
        foldir = dir;
        foldnames = {foldir.name};
        omitlist = regexp(foldnames,'(^\.*|micro|behavioral.*|eeg.*|raw)','match'); % Consider replacing with getDirNamesRegexp.m
        omitlist = [omitlist{:}];
        foldir1 = foldir(~ismember({foldir.name},omitlist));
        foldnames1 = {foldir1.name};
        tocopy1 = strcat(subdir,filesep,sub,filesep,foldnames1); %Add subpath in front of each cell
        sess_unique = unique(sess_list);
        for jj = 1:numel(sess_unique) % Iterate through sessions
            sess = sess_unique{jj};
            barray{subn,1} = sub;
            barray{subn,2} = sess;

            % Get behavioral and session folders to add to copy list
            % (tocopy1)
            for fold2 = 1:numel(foldnames2)
                foldsess = foldnames2{fold2};
                foldpath = [subdir filesep sub filesep foldsess];
                if exist(foldpath,'dir')
                    if contains(foldsess,'behavioral') && jj == 1
                        sess2copy = Task2Copy; % Get the folder named according to the task you want
                    elseif contains(foldsess,'behavioral') && jj > 1 % Not the first session, we already have this in the list
                        continue;
                    elseif contains(foldsess,'eeg')
                        sess2copy = sess;
                    else
                        fprintf('%s\n',['ERROR!! foldsess: ' foldsess ' does not contain behavior or eeg!!!']);
                        keyboard
                    end
                    foldsesspath = [foldpath filesep sess2copy];
                    if exist(foldsesspath,'dir')
                        tocopy1{end+1} = foldsesspath;
                    else
                        fprintf('%s\n',[foldsesspath ' does not exist.']);
                    end

                else
                    fprintf('%s\n',['ERROR!! ' foldpath ' does not exist!!!']);
                    keyboard
                end
            end
            
            
        end
    end

    fprintf('\t%s\n','Directories to copy: ');
    fprintf('\t\t%s\n',tocopy1{:});
    tocopy_all{subn} = tocopy1;
end
fprintf('\n\n%s\n\n','!!!!!!!!!!!!!!!!! CHECK THE ABOVE OUTPUT BEFORE PROCEEDING WITH STEP 3 !!!!!!!!!!!!!!!!!');

%% STEP 3: COPY

outID = fopen([outputdir filesep 'copy_info_' date '.txt'],'w');
fprintf(outID,'%s\n\n\n',['Copy report for ' num2str(numel(sublist)) ' subjects and data copied for Naming Task: ']);

tstart = tic;

for subn2 = 1:numel(sublist)
    sub2 = sublist{subn2};
    outsubdir = [outputdir filesep sub2];
    if ~exist(outsubdir,'dir')
        mkdir(outsubdir)
    end
    fprintf(outID,'%s\n',[sub2 ': ']);
    fprintf(outID,'\t%s\n','Directories to copy: ');
    tocopy_list = tocopy_all{subn2};
    if isempty(tocopy_list)
        fprintf(outID,'\t\t%s\n',['Nothing to copy for Subject ' sub2 '. Skipping to next Subject.']);
        rmdir(outsubdir);
        continue;
    end
        
    for ii = 1:numel(tocopy_list)
        tocopy = tocopy_list{ii};
        outfolddir = regexp(tocopy,'NIH\d{3}/.*(?=/)','match');
        outfoldir = [outfolddir{:}];
        if ~isempty(outfolddir) % Make any new directories
            outfoldpath = [outputdir filesep outfoldir];
            if ~exist(outfoldpath,'dir')
                mkdir(outfoldpath)
            end
        end
        
        toprint = tocopy(regexp(tocopy,'NIH\d.*'):end);
        fprintf(outID,'\t\t%s',toprint);
        fprintf(outID,'\t%s',' ... ');
        
        if exist(tocopy) % Copy
            outcopy = [outputdir filesep toprint];
            if ~exist(outcopy,'dir')
                status = copyfile(tocopy,outcopy);
                if status == 1
                    fprintf(outID,'\t%s\n','Successful.');
                else
                    fprintf(outID,'%s\n','ERROR!!! Tried to copy but unsuccessful!!!');
                end
            else
                fprintf(outID,'%s\n',['ERROR!!! Output directory already exists!!! (' outcopy ')']);
            end         
        else
            fprintf(outID,'%s\n','ERROR!!! Input directory does not exist!!!');
        end
    end
    fprintf(outID,'\n');
end
    
tstop = toc(tstart);

fprintf(outID,'%s\n','Copying Complete.');
fprintf(outID,'Elapsed Time: %d min, %f s\n', floor(tstop/60), rem(tstop,60));

fclose all;

fprintf('%s\n','Done. Please refer to copy_info.txt for copying information.');