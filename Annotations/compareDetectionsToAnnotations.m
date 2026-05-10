function [sens,falseDetections] = compareDetectionsToAnnotations(sl)

timeWindow = 5; % Minutes
timeWindow = timeWindow * 60 * sl.Fs; % Samples

%% Pull detections 

sl.populateSpikes('forceNew',true);

[iiAuto,jjAuto] = find(sl.spikeDetectionResults.rasters); % In samples, from clip start...

%% Pull annotations 

[s, clipDetails] = spreadsheetToSDTimes(sl.rootFolder,sl.subj,'chanNames',sl.chanNames);
populateSeqFromAnnotations(sl,s,clipDetails);

[iiManual,jjManual] = find(sl.spikeDetectionResults.rasters); 

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


end