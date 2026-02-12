function util_set_data_env(enviroType)
% DEPRECATED - DO NOT USE

% Desc: set data environment variables
% Inputs: 
%   enviroType: where this directory lives
% Outputs:
%   environment variables contained within direnv.enviroType.csv
% Example:
%   enviroType = '~/Desktop/ws-cocjin/ROIstore/';
%   util_set_data_env
% Note to self:
%   To do:
%       [home]/params/ 
%%

% if there is no input given, dir list will be pulled from
% direnv.enviroType.csv in opts/environment folder
if nargin == 0
    enviroType = 'default';
end

% grab environment variable matches
dirCsvPath = [getenv('LR_DIR_ENV') '/data.' enviroType '.csv'];

if isempty(dirCsvPath)
    fprintf('Error: No data.%s.csv file found. Follow these steps:\n', enviroType);
    fprintf('\tCreate a folder named environments on your local computer\n');
    fprintf('\tCopy to environments the templates from eeg_toolbox/localize/other/environment_templates.\n');
    fprintf('\tRename these files to remove the template_ prefix\n');
    fprintf('\tFill out data.environments.csv rows to point to the root eeg/ dir you will be working in.\n');
    
    error('Did not find data.%s.csv. Ensure environments folder is on your path.', enviroType);
end
dirList = readtable(dirCsvPath);

% set environment variables
disp('----DATASTORE ENVIRONMENT VARIABLES-----');
for iEnv = 1:height(dirList)
    % get mappings
    myEnv = dirList{iEnv,'env'}{:}; 
    myDir = dirList{iEnv,'dir'}{:};
    
    % set mappings
    setenv(myEnv, myDir);
   disp([myEnv ' = ' getenv(dirList{iEnv,'env'}{:})])
end
disp('----------------------------------------');
return