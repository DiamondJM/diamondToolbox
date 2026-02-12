

function mySingleton = localizationManager(mySingleton)
%% Background basics 

sourceLoc = mySingleton.sourceLoc;
paramStruct = struct;

propagationSpeed = 300;
sensorDistance = 12;
filterWindow = [2 12];

paramStruct.propagationSpeed = propagationSpeed;
paramStruct.sensorDistance = sensorDistance;
paramStruct.filterWindow = filterWindow; 

sourceLoc.paramStruct = paramStruct; 

sourceLoc = getSensors(sourceLoc);

lengthTs = size(sourceLoc.TimeSeries,1);
sourceLoc.lengthTs = lengthTs;


%% Extract phase and power information 

dataMaster = extractPhasePow(sourceLoc);

[dataMaster, params] = postProcessPhase(dataMaster);
paramStruct = loadParams(params,paramStruct);

%% Compute Radius differences over leads 

dPosition = computeDPosition(sourceLoc,dataMaster);

%% Localize 

[localizationResults,params] = localizationFunction(sourceLoc,dPosition);
paramStruct = loadParams(params,paramStruct);


%% Wrap up

sourceLoc.localizationResults = localizationResults; 

sourceLoc.paramStruct = paramStruct;
mySingleton.sourceLoc = sourceLoc;


end
