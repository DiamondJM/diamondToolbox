function pathList = util_list_dirs(varargin)
% pathList = util_list_dirs('rootEEGdir',rootEEGdir [,'path_focus','proc'])
% key-value parameter pairs: 
%  'rootEEGdir', value - pass in the root eeg dir. eg /FRNU/data/eeg/
%  'path_focus', value [default='proc'] - pass in one of these: subj,best,proc,standalone,sink

%% parse inputs
ip = inputParser;
ip.addParameter('path_focus','proc');
ip.addParameter('rootEEGdir','');
ip.parse(varargin{:});
path_focus = ip.Results.path_focus;
rootEEGdir = ip.Results.rootEEGdir;

%%
%filesPath = which('util_dir_tree.csv');
pathListRaw = util_priv_list_datatree(path_focus); % default to subject facing for time being
switch path_focus
    case {'subj', 'best'}
        LRdir = turnary(isempty(rootEEGdir), getenv('LR_DIR_SUBJFOC'), rootEEGdir);
        pathList = structfun(@(x) GetFullPath(fullfile(LRdir,x)),pathListRaw,'UniformOutput',0);
    case 'proc'
        LRdir = turnary(isempty(rootEEGdir), getenv('LR_DIR_PREPROC'), rootEEGdir);
        pathList = structfun(@(x) GetFullPath(fullfile(LRdir,x)),pathListRaw,'UniformOutput',0);
    case 'standalone'
        LRdir = turnary(isempty(rootEEGdir), getenv('LR_DIR_PREPROC'), rootEEGdir);
        pathList = structfun(@(x) GetFullPath(fullfile(LRdir,x)),pathListRaw,'UniformOutput',0);
    case 'sink'
        LRdir = turnary(isempty(rootEEGdir), getenv('LR_DIR_PREPROC'), rootEEGdir);
        pathList = structfun(@(x) GetFullPath(fullfile(LRdir,x)),pathListRaw,'UniformOutput',0); % for now. will eventually switch to its own dir. maybe live on server.            
end
if isempty(LRdir)
    warning('Missing rootEEGdir parameter or environment variable')
end
return
