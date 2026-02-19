function util_set_static_env
% DEPRECATED DO NOT USE

% Desc: set static environment variables
% Inputs: 
%   optsPath: where this directory lives
% Outputs:
%   environment variables contained within LR_ABS.csv
% Example:
%   optsPath = '~/Desktop/ws-cocjin/ROIstore/';
%   util_set_static_env(optsPath)
% Note to self:
%  'util_set_static_env.m' fuck you, greedy ass matlab

%%
% grab environment variable matches
dirCsvPath = which('globals.csv');
if isempty(dirCsvPath)
    fprintf('Error: No globals.csv file found. Follow these steps:\n');
    fprintf('\tCreate a folder named environments on your local computer\n');
    fprintf('\tCopy to environments the templates from eeg_toolbox/localize/other/environment_templates.\n');
    fprintf('\tRename these files to remove the template_ prefix\n');
    fprintf('\tFill out globals:environment and code rows to point to your toolbox and environment folders.\n');
    
    error('Did not find globals.csv. Ensure environments folder is on your path');
end

dirList = readtable(dirCsvPath);

% set environment variables
disp('------STATIC ENVIRONMENT VARIABLES------');
for iEnv = 1:height(dirList)
    % get mappings
    myEnv = dirList{iEnv,'env'}{:}; 
    myDir = dirList{iEnv,'val'}{:};
    
    % set mappings
    setenv(myEnv, myDir);
  disp([myEnv ' = ' getenv(dirList{iEnv,'env'}{:})])
end
disp('----------------------------------------');
return