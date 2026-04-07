function tempScript(allPts)

rootFolder = '/Users/diamondjm/Documents/MATLAB/rootFolder';

for ii = 1:length(allPts)

    thisPt = sprintf('CSSD0%d',allPts(ii)); 

    load(fullfile(rootFolder,thisPt,'sl.mat'),'sl');

    edfFile = dir(fullfile(sl.subjFolder, 'ts', '*.edf'));
    edfPath = fullfile(edfFile(1).folder, edfFile(1).name);


    [timeSeries, Fs, ~] = sl.loadTsFromEdf(edfPath,sl.chanNames);

    sl.timeSeries = timeSeries; 
    sl.Fs = Fs; 

    sl.disposeOfBadSegment; 
    sl.downsampleTs; 

    save(fullfile(sl.subjFolder,'sl.mat'),'sl')
   

end

end