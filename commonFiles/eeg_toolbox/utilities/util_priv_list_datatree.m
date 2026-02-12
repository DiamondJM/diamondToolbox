function dirList = util_priv_list_datatree(datatreeType,varargin)
datatreePath = which('util_dir_tree.csv');
if nargin == 0
    datatreeType = 'subjDir';
end
switch datatreeType
    case 'proc'
        datatreeType = 'procDir';
    case 'subj'
        datatreeType = 'subjDir';
    case 'standalone'
        datatreeType = 'standAlone';
    case 'sink'
        datatreeType = 'sink';
    case 'best'
        datatreeType = 'bestDir';
end

fuqu_table = readtableSafe(datatreePath);

%% construct dirList from table
altName = {};
dirList = struct();
for iRow = 1:height(fuqu_table)
    myEntry = fuqu_table{iRow,datatreeType}{:};
    if strcmp(myEntry,'alt')
            dirList.(fuqu_table{iRow,'id'}{:}) = fuqu_table{iRow,'alt'}{:};
            altName = [altName fuqu_table{iRow,'id'}{:}];
    else
            dirList.(fuqu_table{iRow,'id'}{:}) = myEntry;
    end
    
end

%% format to defined pattern where [base] is dirList.base, [treedir] is dirList.[], and [subjId] is provided eventually
% define pattern
origpattern = dirList.pattern;
mybase = dirList.base;
mypattern = strrep(origpattern, '[base]',mybase);

% remove pattern and base from dirList
dirList = rmfield(dirList, 'pattern');
dirList = rmfield(dirList, 'base');

% format each directory path based on pattern
dirNames = fieldnames(dirList); % augment names
for iDir = 1:length(dirNames)
    myDirName = dirNames{iDir};
    if strcmp(dirList.(myDirName),'x') || ismember(myDirName, altName) 
    else
        dirList.(myDirName) = strrep(mypattern,'[treedir]',dirList.(myDirName));
    end
end

%% get full path, if so desired
if length(varargin) > 0
    if strcmp(varargin{1},'full')
        switch datatreeType
            case 'subjDir'
                dirList = structfun(@(x) fullfile(getenv('LR_DIR_SUBJFOC'),x),dirList,'UniformOutput',0);
            case 'procDir'
                dirList = structfun(@(x) fullfile(getenv('LR_DIR_PREPROC'),x),dirList,'UniformOutput',0);
            case 'standAlone'
                dirList = structfun(@(x) fullfile(getenv('LR_DIR_PREPROC'),x),dirList,'UniformOutput',0);
            case 'sink'
                dirList = structfun(@(x) fullfile(getenv('LR_DIR_PREPROC'),x),dirList,'UniformOutput',0);
        end
    end
end
return