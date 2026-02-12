function rerunAllSpikesSeizures(varargin) 

fprintf('Patient data will be coming from %s.\n',pullServerDirectory); 
pause(1); 


p = inputParser;
addParameter(p,'allPtSz',nan);
addParameter(p,'allPtBl',nan);
addParameter(p,'hardwareType','subdural'); 
parse(p,varargin{:})
allPtSz = p.Results.allPtSz;
allPtBl = p.Results.allPtBl;

isHopeless = {'NIH060','NIH076','NIH085'};

%%

if ismember('allPtSz',p.UsingDefaults)
    szTable = load_SZ_info([],'/Users/diamondjm/Documents/MATLAB/seizureTable/szTableResearch.xlsx');
    allPtSz = unique(szTable.subjID);
    
    allPtSz = setdiff(allPtSz,isHopeless); 
    
end

% szExist = determinePatientsFromDirectory(pullSeizureDirectory);
% allPtSzNum = rectifySubj(allPtSz); 
% [~,inds] = setdiff(allPtSzNum,szExist); 
% allPtSz = allPtSz(inds); 

cd(pullSeizureDirectory); 

if exist('ptResultsSz.mat','file')
    load('ptResultsSz','ptResultsSz');
else
    ptResultsSz = cell(size(allPtSz)); 
end

for ii = 1:length(allPtSz) 
    if isa(ptResultsSz{ii},'char') && contains(ptResultsSz{ii},'Success'); continue; end 
    try 
        eventManager(allPtSz(ii)); 
        ptResultsSz{ii} = sprintf('Success for patient %s.',allPtSz{ii}); 
    catch e
        ptResultsSz{ii} = e; 
        fprintf('Failure for %s.\n',allPtSz{ii}); 
    end
    
    save('ptResultsSz','ptResultsSz'); 
end


%%

if ismember('allPtBl',p.UsingDefaults)
    blTable = load_BL_info([],'/Users/diamondjm/Documents/MATLAB/seizureTable/blTableResearch.xlsx');
    allPtBl = unique(blTable.subjID);
    
    allPtBl = setdiff(allPtBl,isHopeless); 
    
end

cd(pullSpikesDirectory); 
% populateSpikesWrapper(allPtBl)

% if exist('ptResultsBl.mat','file')
%     load('ptResultsBl','ptResultsBl');
% else
%     ptResultsBl = cell(size(allPtSz)); 
% end

ptResultsBl = cell(size(allPtBl));

alreadyDone = determinePatientsFromDirectory('/Users/diamondjm/Documents/Josh_work/sourceLocSpikes/2023Feb25/myStructs/localizedStruct'); 

for ii = 1:length(allPtBl)

    if ismember(rectifySubj(allPtBl{ii}),alreadyDone)
        ptResultsBl{ii} = 'AlreadyDone'; 
        continue
    end
    
    
    try 
        populateSpikesWrapper(allPtBl(ii)); 
        ptResultsBl{ii} = sprintf('Success for patient %s.',allPtBl{ii});
    catch e
        ptResultsBl{ii} = e; 
        fprintf('Failure for %s.\n',allPtBl{ii}); 
    end
    save('ptResultsBl','ptResultsBl'); 
end


end