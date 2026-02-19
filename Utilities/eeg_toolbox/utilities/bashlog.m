function success = bashlog(command, logFile, errorFile, errorKeyword, append)
% BASHLOG executes a bash shell command and redirects output to a log. 
%
% success = bashlog(command, logFile, errorFile, errorKeyword, append)
%
% All output (stdout and stderr) is redirected to logFile, while 
% output lines containing the (case-insensitive) word errorKeyword are also 
% logged to the errorFile and output to console. If there are no errors, 
% the the command returns true. The errorFile contains lines with the error
% keyword (plus the subsequent 0 lines).
%
% errorKeyword may be a string or a cell array of strings.
%
% If append is passed false, logFile and errorFile are moved to 
% trash if they already exist!

EXTRA_LINES = 0; % error line + X following lines

if ~exist('append', 'var') % default is append
    append = true;
end

% get rid of existing files
if append
    appendFlag = '-a ';
    eof = getEOF(errorFile);
else
    appendFlag = '';
    if exist(errorFile, 'file')
        movefile(errorFile, '~/.Trash');
    end
    if exist(logFile, 'file')
        movefile(logFile, '~/.Trash');
    end
end
    
if iscellstr(errorKeyword)
    word = strjoin(errorKeyword,' or ');
else
    word = errorKeyword;
end
% print some explanation to the user
fprintf(['\nRunning a terminal command now. You may monitor progress via '...
    'terminal with command "tail -f [logfile]."\n']);
fprintf('Output containing words "%s" will be considered an error and output to', word);
fprintf(' a log file and to screen.\n');
fprintf('-----------------------------------------------------------------\n');
fprintf(' Command: %s\n', command);
fprintf('Full Log: %s\n', logFile);
fprintf('  Errors: %s\n', errorFile);
fprintf('Running... [started %s]\n', char(datetime('now')));
fprintf('-----------------------------------------------------------------\n');

% build a grep search string for keywords
if iscellstr(errorKeyword)
    grepKeyword = sprintf('grep -iE -A %d "%s"', ...
        EXTRA_LINES, strjoin(strtrim(errorKeyword), '|'));
else
    grepKeyword = sprintf('grep -i -A %d %s', EXTRA_LINES, errorKeyword);
end

% execute the command
pipedCommand = sprintf('%s 2>&1 | tee %s%s | %s | tee %s%s',...
    command, appendFlag, logFile, grepKeyword, appendFlag, errorFile);
system(pipedCommand);

% check if the error file is empty (or nothing was added)
if ~exist(errorFile, 'file') || isBlank(errorFile)
    success = true;
elseif append && exist('eof', 'var')
    success = (eof == getEOF(errorFile)); % no errors added
else
    success = false;
end

if success
    fprintf('Finished successfully! [%s]\n', char(datetime('now')));
else
    fprintf('Finished wih at least one error keyword in output. [%s]\n', char(datetime('now')));
end

end % bashlog function

function isblank = isBlank(filename)
    isblank = 1;
    if ~exist(filename, 'file'), return; end
    fd = fopen(filename, 'r');
    fseek(fd, 0, 'bof');
    bof = ftell(fd);
    fseek(fd, 0, 'eof');
    eof = ftell(fd);
    isblank = (bof == eof);
    fclose(fd);
end

function eof = getEOF(filename)
    eof = 0;
    if ~exist(filename, 'file'), return; end
    fd = fopen(filename, 'r');
    fseek(fd, 0, 'eof');
    eof = ftell(fd);
    fclose(fd);
end