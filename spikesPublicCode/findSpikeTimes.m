function sourceLoc = findSpikeTimes(inputData,peak_win, z_thresh, cts_thresh,varargin)


%% parse inputs
p = inputParser;

% specific to this code
def_amp_scale = 3;
addParameter(p,'amp_scale',def_amp_scale,@(x) isscalar(x));

% For passing in time series. 
addParameter(p,'maxNegPeakWidth',50); 
addParameter(p,'maxPosPeakWidth',Inf); 
addParameter(p,'peakSeparation',0); 
addParameter(p,'trackPeaks',false); 
addParameter(p,'maxPeakHeight',Inf); 

parse(p,varargin{:})
amp_scale=p.Results.amp_scale;

maxNegPeakWidth = p.Results.maxNegPeakWidth;
maxPosPeakWidth = p.Results.maxPosPeakWidth;
trackPeaks = p.Results.trackPeaks; 
peakSeparation = p.Results.peakSeparation; 
maxPeakHeight = p.Results.maxPeakHeight; 

if trackPeaks
    posPeakSeparation = peakSeparation; negPeakSeparation = 0; 
    posPeakHeight = z_thresh; negPeakHeight = 0; 
else 
    negPeakSeparation = peakSeparation; posPeakSeparation = 0; 
    posPeakHeight = 0; negPeakHeight = z_thresh; 
end

%% load data and set up parameters/arrays

if ~isstruct(inputData); error('Check input arguments. X_raw should be a struct now, with Fs field.'); end
Fs = inputData.Fs;
chanNames = inputData.chanNames;
timeSeries = inputData.timeSeries;

dims=size(timeSeries);

hp=floor(peak_win*Fs); % window around peak

%% find peaks based on polarity
warning('off','signal:findpeaks:largeMinPeakHeight');

fullRaster = sparse(dims(1),dims(2));
waveformsMaster = cell(1,dims(2)); 

parfor(kk=1:length(chanNames))
    
    currentRaster = sparse(dims(1),1);
    waveforms = cell(size(currentRaster)); 

    %% Peaks and troughs 
    
    xCurrent = timeSeries(:,kk);
    [pHeight,pInd]=findpeaks(xCurrent,'MinPeakprominence',z_thresh,'MinPeakHeight',posPeakHeight,'MaxPeakWidth',maxPosPeakWidth,'minPeakDistance',posPeakSeparation);

    isBad = pHeight > maxPeakHeight; 
    pHeight(isBad) = []; 
    pInd(isBad) = []; 
    
    [nHeight,nInd]=findpeaks(-xCurrent,'MinPeakProminence',z_thresh,'MinPeakHeight',negPeakHeight,'MaxPeakWidth',maxNegPeakWidth,'minPeakDistance',negPeakSeparation);
    
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
        inds = matchingMatrix >= -hp & matchingMatrix <= hp;
        
        [pFind,nFind] = find(inds);
        
        heightMatrix = false(size(pFind));
        for jj = 1:length(pFind)
            heightMatrix(jj) = pHeight(pFind(jj)) + Nh2Buffer(nFind(jj),ii) >= amp_scale * z_thresh;
        end
        
        nFind = nFind(heightMatrix); 
        pFind = pFind(heightMatrix); 
                
        % Waveform
        if trackPeaks
            for jj = 1:length(pFind)
                ts = nan(1,2 * hp + 1);
                % tsBb = ts; 
                
                startInd = max(pInd(pFind(jj)) - hp,1);
                startBuffer = max(-(pInd(pFind(jj)) - hp) + 2,1);
                
                endInd = min(pInd(pFind(jj)) + hp,dims(1));
                endBuffer = min(hp-(pInd(pFind(jj))-dims(1)) + 1,2 * hp + 1);
                % ts = xCurrent(pInd(pFind(jj))- hp:pInd(pFind(jj)) + hp);
                ts(startBuffer:endBuffer) = xCurrent(startInd:endInd);
                waveforms{pInd(pFind(jj))} = ts;
            end
        else
            for jj = 1:length(nFind)
                
                ts = nan(1,2 * hp + 1);
                % I'll start by filling this with NaNs, for the rare event that
                % we have spikes at the very beginning of the time series.
                % tsBb = ts; 
                
                
                startInd = max(nIndBuffer(nFind(jj),ii) - hp,1);
                startBuffer = max(-(nIndBuffer(nFind(jj),ii) - hp) + 2,1);
                
                endInd = min(nIndBuffer(nFind(jj),ii) + hp,dims(1));
                endBuffer = min(hp-(nIndBuffer(nFind(jj),ii)-dims(1)) + 1,2 * hp + 1);
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

badCounts = sum(fullRaster) < cts_thresh; 
fullRaster(:,badCounts) = false; 

%% Pack up 

sourceLoc.rasters = fullRaster;
sourceLoc.waveforms = waveformsMaster;


end

