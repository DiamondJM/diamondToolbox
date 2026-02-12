function pathList = util_list_paths(varargin)
ip = inputParser;
ip.addParameter('path_focus','proc');
ip.addParameter('rootEEGdir','');
ip.parse(varargin{:});
path_focus = ip.Results.path_focus;
rootEEGdir = ip.Results.rootEEGdir;

% put in if-then statement some time
filesPath = which('util_path_tree.csv');

switch path_focus
    case {'subj', 'best'}
        pathListRaw = pd_loadFilesStatics(filesPath, path_focus);
        LRdir = turnary(isempty(rootEEGdir), getenv('LR_DIR_SUBJFOC'), rootEEGdir);
        pathList = structfun(@(x) fullfile(LRdir,x),pathListRaw,'UniformOutput',0);
    case 'proc'
        pathListRaw = pd_loadFilesStatics(filesPath, 'proc');
        LRdir = turnary(isempty(rootEEGdir), getenv('LR_DIR_PREPROC'), rootEEGdir);
        pathList = structfun(@(x) fullfile(LRdir,x),pathListRaw,'UniformOutput',0);
    case 'standalone'
        pathListRaw = pd_loadFilesStatics(filesPath, 'standalone');
        LRdir = turnary(isempty(rootEEGdir), getenv('LR_DIR_PREPROC'), rootEEGdir);
        pathList = structfun(@(x) fullfile(LRdir,x),pathListRaw,'UniformOutput',0);
    case 'sink'
        pathListRaw = pd_loadFilesStatics(filesPath, 'sink');
        LRdir = turnary(isempty(rootEEGdir), getenv('LR_DIR_PREPROC'), rootEEGdir);
        pathList = structfun(@(x) fullfile(LRdir,x),pathListRaw,'UniformOutput',0); % for now. will eventually switch to its own dir. maybe live on server.
end
return

function pathList = pd_loadFilesStatics(filesoldPath, filesoldType, varargin)
%% construct raw table and file definitions
% load directories and old files definitions
dirList = util_priv_list_datatree(filesoldType); % default to subject facing for time being
pd_table = readtable(filesoldPath);

%% construct pathList for files
pathList= struct();
for iRow = 1:height(pd_table)
    parentdir = dirList.(pd_table{iRow,'dirid'}{:});
    pathList.(pd_table{iRow,'id'}{:}) = GetFullPath(fullfile(parentdir,pd_table{iRow, 'filename'}{:}));
end
%%
return
