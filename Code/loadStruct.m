function [myStruct,failureFlag] = loadStruct(subj,baseDir,varargin)

p = inputParser;
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
useAlternative = p.Results.useAlternative;

if nargin == 1; baseDir = 'myStructs'; end

if useAlternative; baseDir = fullfile(baseDir,'useAlternative'); end 

subj = rectifySubj(subj); 

failureFlag = false;
% baseDir = 'myStructs';
fn = sprintf('myStruct%d.mat',subj);

fprintf('In loading.\n'); 

if exist(fullfile(baseDir,fn),'file') == 2
    load(fullfile(baseDir,fn),'myStruct');
else
    myStruct = pullEventFromServer(subj,'seizure');
    failureFlag = true; 
end 


end