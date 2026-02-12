function [dPosition, params] = computeDPosition(sourceLoc,dataMaster)


%% Preamble 

freqMaster = dataMaster.freq; 
phaseMaster = dataMaster.phase; 
sensors = sourceLoc.sensors.sensorInds;

freqMean = mean(freqMaster,2,'omitnan'); 

% Let's find a bound on frequency difference informed by the doppler
% effect.
propagationSpeed = sourceLoc.paramStruct.propagationSpeed; 
expectedSourceSpeed = 20; 
acceptibleRatio = (propagationSpeed + expectedSourceSpeed) / (propagationSpeed - expectedSourceSpeed);

lengthTs = sourceLoc.lengthTs;
params.acceptibleRatio = acceptibleRatio; 

propagationSpeed = sourceLoc.paramStruct.propagationSpeed;

dPosition = nan(size(sensors,1),lengthTs);

for ii = 1:size(sensors,1)
    currentSensor = sensors(ii,:);
    
    currentPhase = phaseMaster(:,currentSensor);
    currentFreq = freqMaster(:,currentSensor);
    
    validInds = all(~isnan(currentPhase),2);
    
    [s,l] = bounds(currentFreq,2);
    similarFreqLog = l ./ s < acceptibleRatio;
    
    validInds = validInds & similarFreqLog; 
    
    currentPhase = currentPhase(validInds,:); 
    % currentFreq = currentFreq(validInds,:); 
    currentFreq = freqMean(validInds); 
    currentFreq = repmat(currentFreq,1,2); 
    
    
    % currentPhase = unwrap(currentPhase,[],2);
    % Correct, but slower 
    for jj = 1:size(currentPhase,1)
        phaseRange = [-2 * pi + currentPhase(jj,2) currentPhase(jj,2) 2 * pi + currentPhase(jj,2)];
        [~,ind] = min(abs(currentPhase(jj,1) - phaseRange));
        currentPhase(jj,2) = phaseRange(ind);
    end
    
    currentPhase = currentPhase - mean(currentPhase,2);
    
    dPosition(ii,validInds) = -(currentPhase(:,1) * propagationSpeed ./ (2 * pi * currentFreq(:,1)) ...
        - currentPhase(:,2) * propagationSpeed ./ (2 * pi * currentFreq(:,2)));
    
    % Positive distance --> currentSensor(1) is CLOSER than sensor(2)
    
end

fprintf('Finished computing delta Position. \n%.f percent of all computed phase differences used. \n',sum(~isnan(dPosition(:))) / numel(dPosition) * 100);