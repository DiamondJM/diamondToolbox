function comparisonResults = compareDetectionsToAnnotations(subj,rootFolder)

load(fullfile(rootFolder,subj,'sl.mat'),'sl')

timeWindow = 5; % Minutes
timeWindow = timeWindow * 60 * sl.Fs; % Samples

%% Pull detections 

sl.populateSpikes('forceNew',true);

rasterDetected = sl.spikeDetectionResults.rasters; 
[iiAuto,jjAuto] = find(rasterDetected); % In samples, from clip start...

%% Pull annotations 

[s, clipDetails] = spreadsheetToSDTimes(sl.rootFolder,sl.subj,'chanNames',sl.chanNames);
populateSeqFromAnnotations(sl,s,clipDetails);

rasterAnnotated = sl.spikeDetectionResults.rasters; 
[iiManual,jjManual] = find(rasterAnnotated); 

%% Compare

timeMatch = abs(iiManual - iiAuto') <= timeWindow; 
leadMatch = jjManual == jjAuto';

sens = any(timeMatch & leadMatch,2); 
sens = sum(sens) / length(sens); 

falseDetections = ~any(timeMatch & leadMatch); 
falseDetections = sum(falseDetections) / length(falseDetections); 

h = dbstack;
fprintf('[%s] Sensitivity is %.2f%%.\n',h.name,sens * 100)
fprintf('[%s] Ostensibly, %.2f%% of the automatic detections are false.\n',h.name,falseDetections * 100)

comparisonResults = struct( ...
    'rasterDetected',  rasterDetected, ...
    'rasterAnnotated', rasterAnnotated, ...
    'timeWindow',      timeWindow,...
    'sens',sens,...
    'falseDetections',falseDetections);   % samples

sl.plotTimeSeries('spikePlottingMode','fromRaster','comparisonResults',comparisonResults)

end