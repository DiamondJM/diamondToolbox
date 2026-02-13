function sourceLoc = findSpikeTimes(inputData,peakWin, zThresh, ctsThresh,varargin)


%% parse inputs
p = inputParser;

% specific to this code
defAmpScale = 3;
addParameter(p,'amp_scale',defAmpScale,@(x) isscalar(x));

% For passing in time series. 
addParameter(p,'maxNegPeakWidth',50); 
addParameter(p,'maxPosPeakWidth',Inf); 
addParameter(p,'peakSeparation',0); 
addParameter(p,'trackPeaks',false); 
addParameter(p,'maxPeakHeight',Inf); 

parse(p,varargin{:})
ampScale = p.Results.ampScale;

maxNegPeakWidth = p.Results.maxNegPeakWidth;
maxPosPeakWidth = p.Results.maxPosPeakWidth;
trackPeaks = p.Results.trackPeaks; 
peakSeparation = p.Results.peakSeparation; 
maxPeakHeight = p.Results.maxPeakHeight; 

if trackPeaks
    posPeakSeparation = peakSeparation; negPeakSeparation = 0; 
    posPeakHeight = zThresh; negPeakHeight = 0; 
else 
    negPeakSeparation = peakSeparation; posPeakSeparation = 0; 
    posPeakHeight = 0; negPeakHeight = zThresh; 
end

%% load data and set up parameters/arrays

if ~isstruct(inputData); error('Check input arguments. X_raw should be a struct now, with Fs field.'); end
Fs = inputData.Fs;
chanNames = inputData.chanNames;
timeSeries = inputData.timeSeries;

dims=size(timeSeries);

peakWinSamples = floor(peakWin*Fs); % window around peak; samples 

%% find peaks based on polarity
warning('off','signal:findpeaks:largeMinPeakHeight');

fullRaster = sparse(dims(1),dims(2));
waveformsMaster = cell(1,dims(2)); 

parfor(kk=1:length(chanNames))
    
    currentRaster = sparse(dims(1),1);
    waveforms = cell(size(currentRaster)); 

    %% Peaks and troughs 
    
    xCurrent = timeSeries(:,kk);
    [pHeight,pInd]=findpeaks(xCurrent,'MinPeakprominence',zThresh,'MinPeakHeight',posPeakHeight,'MaxPeakWidth',maxPosPeakWidth,'minPeakDistance',posPeakSeparation);

    isBad = pHeight > maxPeakHeight; 
    pHeight(isBad) = []; 
    pInd(isBad) = []; 
    
    [nHeight,nInd]=findpeaks(-xCurrent,'MinPeakProminence',zThresh,'MinPeakHeight',negPeakHeight,'MaxPeakWidth',maxNegPeakWidth,'minPeakDistance',negPeakSeparation);
    
    isBad = false(size(nHeight));
    nHeight(isBad) = [];
    nInd(isBad) = [];


    %% Matching
        
    maxNumel = 1000000; 
    bufferSize = min(100, ceil(maxNumel / length(pInd)));
    
    nIndBuffer = buffer(nInd,bufferSize);    
    Nh2Buffer = buffer(nHeight,bufferSize);
    
    Nh2Buffer(nIndBuffer == 0) = NaN; 
    nIndBuffer(nIndBuffer == 0) = NaN;
    
    for ii = 1:size(nIndBuffer,2)
        
        matchingMatrix = pInd - nIndBuffer(:,ii)';
        
        % inds = matchingMatrix >= -hp & matchingMatrix < 0;
        % The above -- for up-deflection to come before down-deflection.
        inds = matchingMatrix >= -peakWinSamples & matchingMatrix <= peakWinSamples;
        
        [pFind,nFind] = find(inds);
        
        heightMatrix = false(size(pFind));
        for jj = 1:length(pFind)
            heightMatrix(jj) = pHeight(pFind(jj)) + Nh2Buffer(nFind(jj),ii) >= ampScale * zThresh;
        end
        
        nFind = nFind(heightMatrix); 
        pFind = pFind(heightMatrix); 
                
        % Waveform
        if trackPeaks
            for jj = 1:length(pFind)
                ts = nan(1,2 * peakWinSamples + 1);
                % tsBb = ts; 
                
                startInd = max(pInd(pFind(jj)) - peakWinSamples,1);
                startBuffer = max(-(pInd(pFind(jj)) - peakWinSamples) + 2,1);
                
                endInd = min(pInd(pFind(jj)) + peakWinSamples,dims(1));
                endBuffer = min(peakWinSamples-(pInd(pFind(jj))-dims(1)) + 1,2 * peakWinSamples + 1);
                % ts = xCurrent(pInd(pFind(jj))- hp:pInd(pFind(jj)) + hp);
                ts(startBuffer:endBuffer) = xCurrent(startInd:endInd);
                waveforms{pInd(pFind(jj))} = ts;
            end
        else
            for jj = 1:length(nFind)
                
                ts = nan(1,2 * peakWinSamples + 1);
                % I'll start by filling this with NaNs, for the rare event that
                % we have spikes at the very beginning of the time series.
                % tsBb = ts; 
                
                
                startInd = max(nIndBuffer(nFind(jj),ii) - peakWinSamples,1);
                startBuffer = max(-(nIndBuffer(nFind(jj),ii) - peakWinSamples) + 2,1);
                
                endInd = min(nIndBuffer(nFind(jj),ii) + peakWinSamples,dims(1));
                endBuffer = min(peakWinSamples-(nIndBuffer(nFind(jj),ii)-dims(1)) + 1,2 * peakWinSamples + 1);
                % ts = xCurrent(nIndBuffer(nFind(jj),ii)- hp:nIndBuffer(nFind(jj),ii) + hp);
                ts(startBuffer:endBuffer) = xCurrent(startInd:endInd);
                
                %             clf; plot(nIndBuffer(nFind(jj),ii)- hp:nIndBuffer(nFind(jj),ii) + hp,ts); hold on
                %             plot(nIndBuffer(nFind(jj),ii),xCurrent(nIndBuffer(nFind(jj),ii)),'ro')
                %             plot(pInd(pFind(jj)),xCurrent(pInd(pFind(jj))),'bo');
                % pause
                
                waveforms{nIndBuffer(nFind(jj),ii)} = ts;
                                
            end
        end
        
        if trackPeaks; currentRaster(pInd(pFind)) = true; % To retain peaks 
        else; currentRaster(nIndBuffer(nFind,ii)) = true;
        end
        
    end
    waveforms(~currentRaster) = []; 
    waveforms = cell2mat(waveforms); 
    waveformsMaster{kk} = waveforms; 
    
    fullRaster(:,kk) = currentRaster;
    
    
end

%% Narrow by counts 

badCounts = sum(fullRaster) < ctsThresh; 
fullRaster(:,badCounts) = false; 

%% Pack up 

sourceLoc.rasters = fullRaster;
sourceLoc.waveforms = waveformsMaster;


end

