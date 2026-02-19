function status = bashlog2(command, fname_log, append)
% bashlog executes a bash command and tee's the output to log file
%
% If append is given and is false, old log file is first deleted if it already exists.

    if nargin < 2
        fname_log = 'log';
    end
    
    if nargin < 3
        append = true;
    end
    
    if append
        append_flag = '-a';
    else
        append_flag = '';
    end

    % Redirect stderr to stdout, copy stdout to logfile
    % Prefix with date and command string
    full_command = sprintf('echo `date` | tee %s %s; echo Command: "%s" | tee -a %s;', ...
        append_flag, fname_log, command, fname_log);
   
    full_command = strcat(full_command, sprintf('%s 2>&1 | tee -a %s', command, fname_log));
    
    status = unix(full_command);

    fprintf('Output logged to %s\n', fname_log);

end
