function [localizationResults,params] = localizationFunction(sourceLoc,dPosition)

%% Original version 

tic

lengthTs = sourceLoc.lengthTs;
sensors = sourceLoc.sensors.sensorInds;

myBp = sourceLoc.braindata.myBp;

geodesic = sourceLoc.geodesic;
geodesicMaster = geodesic.geodesicDistances;
isLeftInds = geodesic.isLeftInds;
isLeftSensors = isLeftInds(sensors(:,1));

lVertices = myBp.surfaces.pial_lh.vertices;
rVertices = myBp.surfaces.pial_rh.vertices;

marginError = .5;
% This is measured in mm. It's the difference between perceived and actual
% delta distance, among our candidate vertices, which we tolerate. 
% Think of this as the thickness of the tolerance band. 

distancesMaster = sourceLoc.sensors.distancesMaster;

distanceNorm = zeros(length(sensors),1);
ddMaster = zeros(size(sensors,1),size(geodesicMaster,2));
boundsMaster = false(size(sensors,1),size(geodesicMaster,2));

distanceThresh = 30;  % mm
% Beyond this distance, candidate source vertices are disregarded. 
% Literature (Smith, Schevon, Nat Comm 2016, among other work) suggests
% that ~3cm is the ceiling above which neural signals are unlikely to
% travel. 

for ii = 1:size(sensors,1)
    
    currentSensor = sensors(ii,:);
    
    distanceNorm(ii) = distancesMaster(currentSensor(1),currentSensor(2));
    ddMaster(ii,:) = geodesicMaster(currentSensor(1),:) - geodesicMaster(currentSensor(2),:);
    % Here, a positive value (ii,jj) in ddMaster suggests that, for
    % currentSubsensor = sensors(ii,:), vertex jj is closer to
    % currentSubsensor(2). 
    
    
    boundsMaster(ii,:) = max([geodesicMaster(currentSensor(1),:); geodesicMaster(currentSensor(2),:)]) < distanceThresh;
    
    % Max or min? 
    % For each vertex, we have distance to sensor 1 and distance to sensor
    % 2. 
    % Call this d1 and d2. 
    % If max, both d1 and d2 must be below the threshold. 
    % If min, at least one, and possibly both, will be below threshold.
    
    % This can be thought of a distance from source to sensor such that,
    % beyond this distance, computed localization results are likely
    % meaningless. 
    
    % In terms of physiology, min is probably meaningless. If
    % distanceThresh is a hard boundary, we want both sensors to satisfy
    % it. 
end



%% Localization

localizationResults = nan(3,lengthTs);

subsensorLength = 4;

perceivedActualCutoff = 0.9; 

parfor jj = 1:lengthTs
    
    locTemp = nan(3,1);
    
    % currentSubsensor = find(~isnan(drAll(:,jj)));
    currentSubsensor = find(abs(dPosition(:,jj)) < distanceNorm * perceivedActualCutoff);
    
    %%%%
    % Continue if we have a short subsensor.
    if length(currentSubsensor) < subsensorLength; continue; end
    %%%%
    
        
    isLeftSubsensor = isLeftSensors(currentSubsensor);
    useLeft = round(sum(isLeftSubsensor) / length(isLeftSubsensor));
    
    if useLeft; currentSubsensor = currentSubsensor(isLeftSubsensor);
    else; currentSubsensor = currentSubsensor(~isLeftSubsensor); 
    end
    
    if length(currentSubsensor) < subsensorLength; continue; end
    
    drSubsensor = dPosition(currentSubsensor,jj);
    
    
    eligiblePoints = sum(boundsMaster(currentSubsensor,:)) >= 1; % Union of all source spaces
    eligiblePointsIndex = find(eligiblePoints);
    
        
    % Which hemisphere are we on? 
    if isLeftSensors(currentSubsensor(1)); eligiblePoints = lVertices(eligiblePoints,:);
    else; eligiblePoints = rVertices(eligiblePoints,:); 
    end
    
    sourceToVertex = zeros(length(currentSubsensor), size(eligiblePoints,1),'single');
    
    emptyCheck = false(size(currentSubsensor));
    
    for ii = 1:length(currentSubsensor)
        
        currentDr = drSubsensor(ii);
        
        diffDistance = ddMaster(currentSubsensor(ii),:);
        maxDistance = boundsMaster(currentSubsensor(ii),:);
        
        chosenVertices = abs(diffDistance - currentDr) < marginError;
        chosenVertices = chosenVertices & maxDistance;
        
        if ~any(chosenVertices)
            emptyCheck(ii) = true;
            warning('Empty vertex set found for step %d. \n',jj);
            continue; 
        end

        if isLeftSensors(currentSubsensor(1)); chosenVertices = lVertices(chosenVertices,:);
        else; chosenVertices = rVertices(chosenVertices,:);
        end
        
        pointsDistances = pdist2(chosenVertices,eligiblePoints);
        sourceToVertex(ii,:) = min(pointsDistances);
        
    end
    
    if sum(~emptyCheck) < subsensorLength; continue; end
    sourceToVertex(emptyCheck,:) = [];
    
    
    [minVal,ind] = min(mean(sourceToVertex .^ 2));
    % Within a given sample, mean is equivalent to sum. 
    % Across them, mean may make more sense, because this way we don't
    % penalize localizations involving more subsensors. 
    ind = eligiblePointsIndex(ind);

    isLeft = isLeftSensors(currentSubsensor(1));
    if isLeft; ind = ind * -1; end
    
    locTemp(1) = ind; 
    locTemp(2) = minVal;
    locTemp(3) = jj;
    
    localizationResults(:,jj) = locTemp;
    
    
end

fprintf('Localization completed in %.2f seconds. \n',toc);
fprintf('%.f percent of samples gave rise to a localization, for a total of %d localized points. \n',sum(~isnan(localizationResults(2,:))) / lengthTs * 100,sum(~isnan(localizationResults(2,:))));

localizationResults = localizationResults(:, ~isnan(localizationResults(2,:)));

%% Quality control

% quantileThresh = quantile(localizationResults(2,:),.9);
% localizationResults(:,localizationResults(2,:) > quantileThresh) = []; 


%% Pack up


params.subsensorLength = subsensorLength;
params.marginError = marginError;
params.distanceThresh = distanceThresh; 
params.perceivedActualCutoff = perceivedActualCutoff;



end

