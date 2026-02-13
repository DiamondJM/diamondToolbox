classdef sourceLocalizer

    properties

        subj

        chanNames
        Fs

        braindata = struct('myBd',[],'myBp',[]);

        spikeDetectionResults = struct('rasters',[],'waveforms',[],...
            'paramStruct',struct(...
            'seqWin',0.1,... seconds
            'maxNegPeakWidth',0.05,... seconds
            'peakWin', 0.1, ... % seconds
            'ampScale',3,... z 
            'trackPeaks',false,...
            'zThresh',3)); 

        seqResults = struct('seriesAll',[],'timesAll',[]);
        sensors

        sourceLocalizationResults = struct('localizationResults',[],'roiResults',[],...
            'paramStruct', struct(...
            'propagationSpeed',300,... % mm / s
            'sensorDistance',30,... mm 
            'subsensorLength',3)); 


        % These are defaults and can be adjusted.
        % propagationSpeed is in mm/s and refers to speed of wave
        % propagation.

        % sensorDistance is in mm. This is the distance sensors are apart
        % from each other.

        % subsensorLength = 3;
        % How many electrodes at a time (at minimum) should participate in
        % source localization.

    end

    properties (Hidden = true, Transient = true)

        subjFolder

    end

    properties (Transient = true)

        timeSeries
        geodesic

    end

    properties (Hidden = true)
        rootFolder
    end

    methods

        function self = sourceLocalizer(subj,rootFolder,timeSeries,chanNames,Fs)

            % Inputs:
            %   subj -- Subject name, e.g. NIH032, subj12
            %   rootFolder -- root folder name.
            %       Expected folder structure:
            %           rootFolder/
            %               <subj>/           – e.g. NIH032, NIH034
            %                   tal/          – output of the localization process
            %                   <data files>  – time series data can live here too
            %   timeSeries    – n x m array (n samples, m channels)
            %   chanNames     – m x 1 cell array of channel name strings
            %   Fs            – sampling rate (Hz)
            %
            % Note:
            %   You can build a wrapper that takes a subject name, loads the
            %   relevant timeSeries, chanNames, and Fs, and passes them to spikesMain.

            self.subj = subj;
            self.rootFolder = rootFolder;

            self.timeSeries = timeSeries;
            self.chanNames = chanNames;
            self.Fs = Fs;

            %% Pull braindata

            self = self.retrieveBraindata;

        end

        function subjFolder = get.subjFolder(self)

            subjFolder = fullfile(self.rootFolder,self.subj);

        end


        function self = retrieveBraindata(self,varargin)

            %% Preamble


            p = inputParser;
            addParameter(p,'forceNew',false);
            addParameter(p,'saving',true);
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;
            saving = p.Results.saving;

            if ~forceNew && ~isempty(self.braindata.myBd); return; end

            %%

            assert(exist('braindata2', 'file') == 2, 'braindata2 not found on path. Please check paths. Navigate to diamondToolbox and use addpath(genpath(pwd)).');
            assert(exist('brainplotter', 'file') == 2, 'brainplotter not found on path. Please check paths. Navigate to diamondToolbox and use addpath(genpath(pwd)).');

            subjDir = fullfile(self.subjFolder,'braindata');

            fName = fullfile(subjDir,'bdBpFull.mat');

            if exist(fName','file') == 2 && ~forceNew
                S = load(fName,'myBd','myBp');
                myBp = S.myBp;
                myBd = S.myBd;
            else

                if forceNew && saving; rmdir(subjDir,'s'); end
                % This can be handled by determining whether or not forceNew is passed
                % into localizer_rois.

                fprintf('Creating braindata objects for %s.\n',self.subj);
                % root = '/Volumes/Shares/FRNU/dataWorking/sz';

                myBd = braindata2(self.subj,self.rootFolder);
                myBp = ez_get_plotter(myBd);

                if saving
                    mkdir(fullfile(subjDir));
                    save(fullfile(subjDir,'bdBpFull'),'myBd','myBp')
                end


            end

            self.braindata.myBp = myBp;
            self.braindata.myBd = myBd;


        end

        function self = populateSpikes(self,varargin)


            p = inputParser;
            addParameter(p,'forceNew',false);
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;

            if ~isempty(self.spikeDetectionResults.rasters) && ~forceNew; return; end

            %%

            fprintf('Calling spike detector with the following parameters. These are adjustable.\n'); 
            disp(self.spikeDetectionResults.paramStruct); 

            % peakWin = 0.1;
            % zThresh = 3;
            % ampScale = 3;
            % maxNegPeakWidth = 50; % ms
            % maxNegPeakWidth = maxNegPeakWidth / 1000 * self.Fs; % Samples
            % trackPeaks = false;

            self = self.findSpikeTimes;


        end

        function self = findSpikeTimes(self)


            %% parse inputs
            % p = inputParser;

            % specific to this code
            % defAmpScale = 3;
            % addParameter(p,'ampScale',defAmpScale,@(x) isscalar(x));

            % For passing in time series.
            % addParameter(p,'maxNegPeakWidth',50);
            % addParameter(p,'maxPosPeakWidth',Inf);
            % addParameter(p,'peakSeparation',0);
            % addParameter(p,'trackPeaks',false);
            % addParameter(p,'maxPeakHeight',Inf);
            % 
            % parse(p,varargin{:})
            % ampScale = p.Results.ampScale;
            % 
            % maxNegPeakWidth = p.Results.maxNegPeakWidth;
            % maxPosPeakWidth = p.Results.maxPosPeakWidth;
            % trackPeaks = p.Results.trackPeaks;
            % peakSeparation = p.Results.peakSeparation;
            % maxPeakHeight = p.Results.maxPeakHeight;

            % Kill the input parser 

            maxNegPeakWidth = self.spikeDetectionResults.paramStruct.maxNegPeakWidth; 
            maxPosPeakWidth = Inf; 
            trackPeaks = self.spikeDetectionResults.paramStruct.trackPeaks;
            ampScale = self.spikeDetectionResults.paramStruct.ampScale; 
            peakSeparation = 0; 
            maxPeakHeight = Inf; 
            ctsThresh = 0; 
            zThresh = self.spikeDetectionResults.paramStruct.zThresh;
            peakWin = self.spikeDetectionResults.paramStruct.peakWin;

            % If you choose to mess around with ctsThresh or other 
            % 'non-reported' parameters, feel free to pass them back with
            % self.spikeDetectionResults.paramStruct.

            if trackPeaks
                posPeakSeparation = peakSeparation; negPeakSeparation = 0;
                posPeakHeight = zThresh; negPeakHeight = 0;
            else
                negPeakSeparation = peakSeparation; posPeakSeparation = 0;
                posPeakHeight = 0; negPeakHeight = zThresh;
            end

            %% load data and set up parameters/arrays

            thisTs = self.timeSeries; 
            tsDims = size(thisTs);

            peakWinSamples = floor(peakWin*self.Fs); % window around peak; samples
            maxNegPeakWidth = maxNegPeakWidth * self.Fs; % Samples 

            %% find peaks based on polarity

            warning('off','signal:findpeaks:largeMinPeakHeight');

            fullRaster = sparse(tsDims(1),tsDims(2));
            waveformsMaster = cell(1,tsDims(2));

            parfor(kk=1:length(self.chanNames))

                currentRaster = sparse(tsDims(1),1);
                waveforms = cell(size(currentRaster));

                %% Peaks and troughs

                xCurrent = thisTs(:,kk);
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

                            endInd = min(pInd(pFind(jj)) + peakWinSamples,tsDims(1));
                            endBuffer = min(peakWinSamples-(pInd(pFind(jj))-tsDims(1)) + 1,2 * peakWinSamples + 1);
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

                            endInd = min(nIndBuffer(nFind(jj),ii) + peakWinSamples,tsDims(1));
                            endBuffer = min(peakWinSamples-(nIndBuffer(nFind(jj),ii)-tsDims(1)) + 1,2 * peakWinSamples + 1);
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

            %% Volume conduction 

            [fullRaster,waveformsMaster] = self.removeVolCond_fromRaster(fullRaster,waveformsMaster);

            %% Pack up

            self.spikeDetectionResults.rasters = fullRaster;
            self.spikeDetectionResults.waveforms = waveformsMaster;


        end

        function self = computeSequences(self,varargin)

            % Formerly computeSequences_includeDuplicate

            %% Preamble

            p = inputParser;
            addParameter(p,'seriesLength',3);
            addParameter(p,'forceNew',false);
            parse(p,varargin{:})
            seriesLength = p.Results.seriesLength;
            forceNew = p.Results.forceNew;

            if ~forceNew && ~isempty(self.seqResults.seriesAll); return; end

            % seqWin = 0.1;
            % I never really gave the above much thought.
            % But 0.1 is a pretty good number here.
            % We wouldn't really be expecting to localize spikes with sensors greater
            % than 0.3 cm apart, because we wouldn't expect signal to necessarily
            % travel that far anyway.
            % So given a sensor spacing of 3 cm, the maximal interval of time receipt
            % is distance (mm) / propagationSpeed (mm/s)
            % = 30 / 300 mm / mm * s
            % = 0.1 s.


            seqWinSamples = ceil(self.spikeDetectionResults.paramStruct.seqWin * self.Fs);  % Samples
            overlapSamples = floor(seqWinSamples / 2);

            rasters = self.spikeDetectionResults.rasters;

            if isempty(rasters); return; end % Bad session

            totalWindows = buffer(1:size(rasters,1),seqWinSamples,overlapSamples,'nodelay');

            %% On to the function

            [winSpkCnt,~]=buffer(full(sum(rasters,2)),seqWinSamples,overlapSamples,'nodelay');
            winSpkCnt=sum(winSpkCnt);
            hasSeries = find(winSpkCnt >= seriesLength);

            % If nothing


            seriesAll = cell(max(winSpkCnt),length(hasSeries));
            seriesAll(:) = {''};
            timesAll = nan(max(winSpkCnt),length(hasSeries));

            if isempty(hasSeries); return; end

            % For parallel
            timesHolder = nan(max(winSpkCnt),1);
            leadsHolder = cell(max(winSpkCnt),1);
            leadsHolder(:) = {''};

            for jj = 1:length(hasSeries)
                currentWindow = rasters(totalWindows(:,hasSeries(jj)),:);
                [row,col] = find(currentWindow); % row times, col leads
                % currentTimes = totalWindows(row,hasSeries(jj));
                % currentLeads = chanNames(col);

                % Here we can deal with handling of duplicates.
                % That is, an instance when a single electrode has multiple IEDs in a
                % single window.
                % We can use a method more delicate than that used in computeSequences.
                % Essentially under our fundamental assumptions, it doesn't make sense
                % to have duplicates. We assume single impulses are sent from the
                % source. For there to be a duplicate, you would need consecutive
                % barrages with low latency (unlikely) or for the signal to change
                % direction and double back (impossible).
                % The alternative is that there's a single discharge but with a complex
                % waveform, which registers as multiple spikes. The hope is that seqWin
                % is SHORT enough that these end up being called their own sequence.

                % There is a balance where if seqWin is too short, then we'll miss
                % sequences (particularly ones where the wave may travel slowly) which
                % take a long time to traverse the array.
                % But if it's too long, we actually lose VERY LITTLE, if waveforms are
                % simple. If waveforms are complex, we're more likely to incorporate a
                % double spike in the same sequence, where we would have done better to
                % call that a separate sequence.

                % It is possible there are some "double spikes" which are impossible to
                % resolve, if the latency between the two spikes is less than the time
                % it takes for the wave to cross the array. There's not much we can do
                % about this (although there's filtering).
                % Alternatively, sometimes what might register as a "double spike" is
                % actually noise and we don't want that. In other words we'd rather
                % that the noise didn't exist and the double spike were registered as a
                % single.

                % Generally, we care more about making sure the entire sequence is
                % recovered than we do about separating double spikes, so our seqWin
                % will err towards the long side.

                % Regarding double spikes that we do capture, they may be 1) noise, or
                % 2) indicative of a complex waveform. Regardless, we can't afford to
                % make our seqWin shorter, and even for the shortest seqWin, some
                % double spikes will register.

                % Therefore we'll assume these double spikes are noise, and address
                % them by taking the median time.

                % I believe it won't work to take a true median (with rounding). The
                % reason is that the spike times act as a "key" to recover the
                % waveforms in waveformsFromSequence. I get that this system is clumsy
                % and possibly merits revision, but it's what we have now.
                % By virtue of that, let's take a "median" which nonetheless mandates
                % that we pick an element.
                % In fact maybe taking a value, that is, that "closest" to the median,
                % is best, in the event that it's not actually noise.
                % The below method picks arbitrarily (not randomely) which is good. In
                % cases of even ties, it should pick the FIRST before the median.

                [a,b,c] = unique(col);
                if length(a) < seriesLength; continue; end
                currentLeads = self.chanNames(a);
                currentTimes = nan(size(a));
                for ii = 1:length(a)
                    % currentTimes(ii) = median(totalWindows(row(c==ii),hasSeries(jj)));

                    % Force to choose a value
                    vals = totalWindows(row(c==ii),hasSeries(jj));
                    if length(vals) > 1; fprintf('Warning: duplicate time value found.\n'); end
                    % This shouldn't happen in light of the peakSeparation value being
                    % set.
                    [~, ind] = min(abs(vals-median(vals)));
                    currentTimes(ii) = vals(ind);
                end


                [currentTimes,inds] = sort(currentTimes,'ascend');
                currentLeads = currentLeads(inds);
                currentTimes(2:end) = diff(currentTimes);

                % temp1(1:length(currentLeads),jj) = currentLeads;
                % temp2(1:length(currentTimes),jj) = currentTimes;

                currentLeadsHolder = leadsHolder;
                currentLeadsHolder(1:length(currentLeads)) = currentLeads;

                seriesAll(:,jj) = currentLeadsHolder;

                currentTimesHolder = timesHolder;
                currentTimesHolder(1:length(currentTimes)) = currentTimes;

                timesAll(:,jj) = currentTimesHolder;

            end

            seriesAll(:,all(isnan(timesAll))) = [];
            timesAll(:,all(isnan(timesAll))) = [];
            % For sequences that became too short after accounting for duplicate
            % electrodes

            [seriesAll,timesAll] = self.removeDuplicates(seriesAll,timesAll);
            % For duplicate SEQUENCES, big difference

            self.seqResults.seriesAll = seriesAll;
            self.seqResults.timesAll = timesAll;

        end


        function self = localizationManagerSpikes(self,varargin)

            p = inputParser;
            addParameter(p,'plotting',true);
            parse(p,varargin{:})
            plotting = p.Results.plotting;


            %% Preamble


            % Populate spikes
            self = self.populateSpikes;

            % Compute sequences
            self = self.computeSequences;

            % propagationSpeed = 300;
            % sensorDistance = 30;

            % Here, sensorDistance can basically be as far apart as we'd like, barring
            % concerns for the distance at which the signal can travel.
            % We do NOT have aliasing concerns like we did previously.
            % Before, a large sensorDistance meant that we required very low
            % frequencies, otherwise we'd have spatial aliasing.
            % Here, interictal sourceLoc are rare (rather than periodic) events. We assume
            % they arrive rarely, and so aliasing is not an issue. Since frequency of
            % these events are taken to be ~0, there's no upper limit on inter-sensor
            % distance (as dictated by frequency).

            self = self.getSensors;

            %% Collect series

            self = self.buildSeriesTable;

            %% Onto localization

            self = self.localizationFunctionSpikes;

            self = self.locDataToRoi;

            %% Plot

            if ~plotting; return; end
            self.plotSurfFun;

        end

        function self = getSensors(self)

            self = self.loadGeodesic;
            sensorDistance = self.sourceLocalizationResults.paramStruct.sensorDistance;

            geodesicDistancesMaster = self.geodesic.geodesicDistances;
            vertexNums = self.geodesic.vertexNums;

            isLeft = self.geodesic.isLeftInds;

            % distancesMaster = geodesicMaster(:,vertexNums);
            % Unfortunately the above doesn't work anymore, now with Depths.
            % Some elements of vertexNums are NaNs, corresponding to depth electrodes I
            % couldn't localize.
            distancesMaster = nan(length(vertexNums));
            isGood = ~isnan(vertexNums);
            distancesMaster(:,isGood) = geodesicDistancesMaster(:,vertexNums(isGood));

            oppositeSideMask = xor(isLeft,isLeft');
            distancesMaster(oppositeSideMask) = nan;
            % distancesMaster(oppositeSideMask) = Inf;
            % warning('Contralateral sensors Inf.');

            [ii,jj] = find(distancesMaster <= sensorDistance);
            sensorsMaster = [ii jj];

            sensorsMaster = sensorsMaster(ii < jj,:);



            sensorStruct.sensorInds = sensorsMaster;
            sensorStruct.distancesMaster = distancesMaster;
            sensorStruct.sensorDistance = sensorDistance; % Can still save this on the way out

            self.sensors = sensorStruct;


        end

        function self = loadGeodesic(self,varargin)

            %% Preamble

            p = inputParser;
            addParameter(p,'unModified',false);
            addParameter(p,'forceNew',false);
            addParameter(p,'saving',true);
            parse(p,varargin{:})
            unModified = p.Results.unModified;
            forceNew = p.Results.forceNew;
            saving = p.Results.saving;

            %% Load or create, unpack

            if ~isempty(self.geodesic) && ~forceNew
                assert(isequal(length(self.chanNames),length(self.geodesic.leadNames)));
                return;
            end



            geodesicMaster = self.collectGeodesicDistances_master('saving',saving);


            %% Modify

            if unModified
                self.geodesic = geodesicMaster;
                return;
            end

            %% Unpack


            [chanLog,chanInds] = ismember(self.chanNames,geodesicMaster.leadNames);
            % missingLeads = setdiff(chanNames,leadsFromGeodesic);
            missingLeads = self.chanNames(~chanLog);

            if ~all(chanLog)
                warning('Channels submitted that are missing from our records.');

                disp(missingLeads);
                % JD 5/6/2023:
                % The new (improved) paradigm is that geodesicDistances should be
                % created ONLY using leads.csv. In other words, functionality which
                % cross-checks leads with time series (typically Jacksheet leads) with
                % leads with localization, has been removed. We should use only and all
                % the leads in leads.csv to create geodesic. Then pull relevant leads
                % from that 'master' list later.

                % So as of that time, there's no solid reason to programatically
                % re-make geodesic.

            end

            % The above doesn't work anymore because we need to permit tolerance of
            % missing values.
            % Under this new paradigm (3/6/23) we will proceed with localization even
            % if there are channels in jacksheet that failed to localize.

            leadsFromGeodesic = repmat({''},size(self.chanNames));
            leadsFromGeodesic(chanLog) = geodesicMaster.leadNames(chanInds(chanLog));
            % Should match chanNames, minus the missing channels

            geodesicDistances = nan(length(self.chanNames),size(geodesicMaster.geodesicDistances,2));
            geodesicDistances(chanLog,:) = geodesicMaster.geodesicDistances(chanInds(chanLog),:);

            vertexNums = nan(1,length(self.chanNames));
            vertexNums(chanLog) = geodesicMaster.vertexNums(chanInds(chanLog));

            isLeftInds = false(size(self.chanNames));
            isLeftInds(chanLog) = geodesicMaster.isLeftInds(chanInds(chanLog));
            % A bit risque to initialize this is as false. It should be fine, though,
            % because all other fields are nan so we won't be able to do anything with
            % the unlocalized ones (chanNames(~chanLog)).

            [leadLocations,leadLocationsPial] = deal(nan(length(chanInds),3));
            leadLocations(chanLog,:) = geodesicMaster.leadLocations(chanInds(chanLog),:);
            leadLocationsPial(chanLog,:) = geodesicMaster.leadLocationsPial(chanInds(chanLog),:);

            distanceFull = nan(1,length(self.chanNames));
            distanceFull(chanLog) = geodesicMaster.distanceFull(chanInds(chanLog));

            %% Pack

            self.geodesic.geodesicDistances = geodesicDistances;
            self.geodesic.vertexNums = vertexNums;
            self.geodesic.isLeftInds = isLeftInds;
            self.geodesic.leadLocations = leadLocations;
            self.geodesic.leadLocationsPial = leadLocationsPial;
            self.geodesic.leadNames = leadsFromGeodesic;

            self.geodesic.distanceFull = distanceFull;

        end

        function geodesicMaster = collectGeodesicDistances_master(self,varargin)

            %% Preamble

            p = inputParser;
            addParameter(p,'saving',true);
            addParameter(p,'forceNew',false);
            parse(p,varargin{:})
            saving = p.Results.saving;
            forceNew = p.Results.forceNew;

            fn = fullfile(self.subjFolder,sprintf('geodesic_%s.mat',self.subj));

            if ~forceNew && exist(fn,'file')
                load(fn,'geodesicMaster');
                return
            end

            tic

            self = self.retrieveBraindata;
            myBd = self.braindata.myBd; myBp = self.braindata.myBp;

            fprintf('Building geodesic distances for %s. \n',self.subj);

            %% Go on to build the geodesic distances

            surfL = myBp.surfaces.pial_lh;
            surfR = myBp.surfaces.pial_rh;
            numVert = myBd.stdNumVertices;

            leadsDir = fullfile(self.subjFolder,'/tal/leads.csv');
            assert(exist(leadsDir,'file'));
            t = readtable(leadsDir);

            isLeftInds = t.x < 0;


            leadLocations = table2array(t(:,{'x','y','z'}));
            leadNames = t.chanName;
            numLeads = length(leadNames);

            %%%%%%%%%%%%%%
            % Which hemi?
            %%%%%%%%%%%%%%

            distanceThresh = 5; % mm

            leadToVertex = nan(1,numLeads);
            distanceFull = nan(1,numLeads);

            % Left
            leadDistance = pdist2(surfL.vertices,leadLocations(isLeftInds,:));
            [valsLeft,closestVert] = min(leadDistance);
            closestVert(valsLeft >= distanceThresh) = nan;
            leadToVertex(isLeftInds) = closestVert;
            distanceFull(isLeftInds) = valsLeft;

            % Right
            leadDistance = pdist2(surfR.vertices,leadLocations(~isLeftInds,:));
            [valsRight,closestVert] = min(leadDistance);
            closestVert(valsRight >= distanceThresh) = nan;
            leadToVertex(~isLeftInds) = closestVert;
            distanceFull(~isLeftInds) = valsRight;

            leadLocationsPial = nan(numLeads,3);

            isLeftInds = revertToVector(isLeftInds);
            assert(size(isLeftInds,1)==1 && size(leadToVertex,1)==1);

            leadLocationsPial(isLeftInds & ~isnan(leadToVertex),:) = surfL.vertices(leadToVertex(isLeftInds & ~isnan(leadToVertex)),:);
            leadLocationsPial(~isLeftInds & ~isnan(leadToVertex),:) = surfR.vertices(leadToVertex(~isLeftInds & ~isnan(leadToVertex)),:);
            % If use alternative is on, this maps STANDARD leadToVertex, using the ALT
            % surface, to Euclidean space

            %%

            geodesicDistanceMaster = nan(numLeads,numVert);

            parfor ii = 1:numLeads

                if isnan(leadToVertex(ii)); continue; end

                if isLeftInds(ii); geodesicDistanceMaster(ii,:) = myBp.dist_geodesic(surfL, leadToVertex(ii));
                else; geodesicDistanceMaster(ii,:) = myBp.dist_geodesic(surfR, leadToVertex(ii));
                end

            end

            geodesicMaster.geodesicDistances = geodesicDistanceMaster;
            geodesicMaster.vertexNums = leadToVertex;
            geodesicMaster.isLeftInds = isLeftInds;
            geodesicMaster.leadLocationsPial = leadLocationsPial;
            geodesicMaster.leadLocations = leadLocations;  % Dural
            geodesicMaster.leadNames = leadNames;
            geodesicMaster.distanceThresh = distanceThresh;
            geodesicMaster.distanceFull = distanceFull;

            fprintf('Geodesic distances built in %.3f seconds \n',toc)

            if saving; save(fn,'geodesicMaster'); end


        end

        function self = buildSeriesTable(self)

            %% Preamble

            intervalMax = self.sourceLocalizationResults.paramStruct.sensorDistance / self.sourceLocalizationResults.paramStruct.propagationSpeed; % Seconds
            intervalMax = intervalMax * self.Fs;

            % subsensorLength = 3;

            % pairHits = zeros(length(chanNames));

            spkLeads = self.seqResults.seriesAll;
            spkTimes = self.seqResults.timesAll;

            self = self.getSensors;


            %% Prepare ingredients

            % sourceLoc = getSensors(sourceLoc);
            % This should be passed in.

            sensorInds = self.sensors.sensorInds;
            sensorFlip = [sensorInds; fliplr(sensorInds)];

            %% Build table

            seriesTable = nan(size(sensorInds,1),size(spkTimes,2));
            if isempty(seriesTable); return; end

            for jj = 1:size(spkTimes,2)

                currentSeries = spkLeads(:,jj);
                currentTimes = spkTimes(:,jj);
                currentSeries = currentSeries(~isnan(currentTimes));
                currentTimes = currentTimes(~isnan(currentTimes));

                if length(currentTimes) < self.sourceLocalizationResults.paramStruct.subsensorLength; continue; end

                currentTimes(1) = 0;
                for ii = 2:length(currentTimes)
                    currentTimes(ii) = sum(currentTimes(ii - 1:ii));
                end


                [~,leadIndices] = ismember(currentSeries,self.chanNames);
                pairIndices = nchoosek(1:length(currentSeries),2);

                leadPairs = leadIndices(pairIndices);

                [a,b] = ismember(leadPairs,sensorFlip,'rows');
                if sum(a) < self.sourceLocalizationResults.paramStruct.subsensorLength; continue; end

                pairIndices(b > size(sensorInds,1),:) = fliplr(pairIndices(b > size(sensorInds,1),:));
                pairIndices = pairIndices(a,:);

                %     leadPairs = leadIndices(pairIndices);
                %     for ii = 1:size(leadPairs,1)
                %         pairHits(leadPairs(ii,1),leadPairs(ii,2)) = pairHits(leadPairs(ii,1),leadPairs(ii,2)) + 1;
                %     end


                b(b > size(sensorInds,1)) = b(b > size(sensorInds,1)) - size(sensorInds,1);
                b = b(a);

                pairTimes = currentTimes(pairIndices);
                diffTimes = -diff(pairTimes,1,2);
                belowInterval = abs(diffTimes) <= intervalMax;
                if sum(belowInterval) < self.sourceLocalizationResults.paramStruct.subsensorLength; continue; end
                % It will be very unusual to get here, since the <100 ms requirement is
                % enforced by the sequence collection scheme.
                % Unless of course we're using a different propagation speed or
                % something like that

                diffTimes = diffTimes(belowInterval);
                b = b(belowInterval);


                seriesTable(b,jj) = diffTimes; % Samples
                % Lead 1 time - lead 2 time

            end

            % spkTimes(:,all(isnan(seriesTable))) = [];
            % spkLeads(:,all(isnan(seriesTable))) = [];
            % % I used to pass these back as outputs. Not sure if they're needed.
            % % If they are, I should get a "not enough output arguments" error.
            % % Again though, if they are, better to pass them back in params
            % % rather than as outputs.
            %
            % structKey(all(isnan(seriesTable))) = [];

            seriesTable(:,all(isnan(seriesTable))) = [];


            %% Pack up

            % params.subsensorLength = subsensorLength;
            % params.structKey = structKey;
            % params.spikeTimes = spkTimes; % For various reasons I find it's helpful to pass this back in params.
            % params.spikeLeads = spkLeads;

            self.sensors.seriesTable = seriesTable;

            % I think it's good to put this into the sensors struct,
            % because seriesTable should have the same height as
            % sensorInds.

        end

        function self = localizationFunctionSpikes(self,varargin)

            p = inputParser;
            addParameter(p,'forceNew',false);
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;

            %% Preamble

            if ~forceNew && ~isempty(self.sourceLocalizationResults.localizationResults); return; end

            %% Onto the function

            tic

            dPosition = self.sensors.seriesTable;

            dPosition = dPosition / self.Fs; % Seconds
            dPosition = dPosition * self.sourceLocalizationResults.paramStruct.propagationSpeed; % mm

            lengthData = size(dPosition,2);
            sensorInds = self.sensors.sensorInds;

            myBp = self.braindata.myBp;

            self = self.loadGeodesic;

            geodesicDistancesMaster = self.geodesic.geodesicDistances;

            isLeftInds = self.geodesic.isLeftInds;
            isLeftSensors = isLeftInds(sensorInds(:,1));

            lVertices = myBp.surfaces.pial_lh.vertices;
            rVertices = myBp.surfaces.pial_rh.vertices;

            marginError = .5;
            % This is measured in mm. It's the difference between perceived and actual
            % delta distance, among our candidate vertices, which we tolerate.
            % Think of this as the thickness of the tolerance band.

            distancesMaster = self.sensors.distancesMaster;

            distanceNorm = zeros(length(sensorInds),1);
            % ddMaster = zeros(size(sensors,1),size(geodesicMaster,2));
            boundsMaster = false(size(sensorInds,1),size(geodesicDistancesMaster,2));

            distanceThresh = 30;  % mm
            % Beyond this distance, candidate source vertices are disregarded.
            % Literature (Smith, Schevon, Nat Comm 2016, among other work) suggests
            % that ~3cm is the ceiling above which neural signals are unlikely to
            % travel.

            for ii = 1:size(sensorInds,1)

                currentSensor = sensorInds(ii,:);

                distanceNorm(ii) = distancesMaster(currentSensor(1),currentSensor(2));

                % ddMaster(ii,:) = geodesicMaster(currentSensor(1),:) - geodesicMaster(currentSensor(2),:);
                % Perhaps move away from initializing this DD master. It creates a
                % massive array which is then broadcast by the parallel pool.
                % We can define these distances anew for each sensor pair.

                boundsMaster(ii,:) = max([geodesicDistancesMaster(currentSensor(1),:); geodesicDistancesMaster(currentSensor(2),:)]) < distanceThresh;

            end



            %% Localization

            localizationResults = nan(3,lengthData);

            subsensorLength = self.sourceLocalizationResults.paramStruct.subsensorLength;
            perceivedActualCutoff = 0.9;

            parfor jj = 1:lengthData

                locTemp = nan(3,1);

                % currentSubsensor = find(~isnan(drAll(:,jj)));
                currentSubsensor = find(abs(dPosition(:,jj)) < distanceNorm * perceivedActualCutoff);
                % Perceived actual cutoff doesn't affect our source localization procedure itself. Rather, it
                % LIMITS the number of electrode pairs included in our procedure. The
                % reason is that, if dPosition(p,jj) is too close to distanceNorm(p),
                % for some sensor pair p, then the hyperbola starts to look like a
                % hairpin, that is, very eccentrically-bent towards the electrode. This
                % limits the number of available cortical surface points, which can
                % cause the hyperbola to become spotty and could mislead results.

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


                % eligiblePoints = logical(sum(boundsMaster(currentSubsensor,:)));
                eligiblePoints = sum(boundsMaster(currentSubsensor,:)) >= 1; % Union of all source spaces
                % This enforces the requirement that the chosen point at least
                % distanceThresh from at least ONE electrode pair.
                eligiblePointsIndex = find(eligiblePoints);


                % Which hemisphere are we on?
                if isLeftSensors(currentSubsensor(1)); eligiblePoints = lVertices(eligiblePoints,:);
                else; eligiblePoints = rVertices(eligiblePoints,:);
                end

                sourceToVertex = zeros(length(currentSubsensor), size(eligiblePoints,1),'single');

                emptyCheck = false(size(currentSubsensor));

                for ii = 1:length(currentSubsensor)

                    currentDr = drSubsensor(ii);

                    currentSensor = sensorInds(currentSubsensor(ii),:);

                    % diffDistance = ddMaster(currentSubsensor(ii),:);
                    diffDistance = geodesicDistancesMaster(currentSensor(1),:) - geodesicDistancesMaster(currentSensor(2),:);
                    maxDistance = boundsMaster(currentSubsensor(ii),:);
                    % maxDistance = max([geodesicMaster(currentSensor(1),:); geodesicMaster(currentSensor(2),:)]) < distanceThresh;

                    chosenVertices = abs(diffDistance - currentDr) < marginError;
                    chosenVertices = chosenVertices & maxDistance;
                    % Refine the hyperbola so that it extends no more than distanceThresh from
                    % either electrode

                    if ~any(chosenVertices)
                        emptyCheck(ii) = true;
                        warning('Empty vertex set found for step %d. \n',jj);
                        continue;
                    end

                    if isLeftSensors(currentSubsensor(1)); chosenVertices = lVertices(chosenVertices,:);
                    else; chosenVertices = rVertices(chosenVertices,:);
                    end

                    % dT = delaunayn(double(chosenVertices));
                    % [~,sourceToVertex(ii,:)] = dsearchn(double(chosenVertices),dT,double(eligiblePoints));
                    % Slower

                    pointsDistances = pdist2(chosenVertices,eligiblePoints);
                    % For each eligible vertex, find the distance from that point to
                    % all points on the hyperbola

                    sourceToVertex(ii,:) = min(pointsDistances);
                    % For each eligible vertex, find the distance from that point to
                    % the CLOSEST point on the hyperbola

                end

                if sum(~emptyCheck) < subsensorLength; continue; end
                sourceToVertex(emptyCheck,:) = [];

                [minVal,ind] = min(mean(sourceToVertex .^ 2)); % L2 norm
                % [minVal,ind] = min(mean(sourceToVertex)); % L1 norm

                ind = eligiblePointsIndex(ind);

                isLeft = isLeftSensors(currentSubsensor(1));
                if isLeft; ind = ind * -1; end


                locTemp(1) = ind;
                locTemp(2) = minVal;
                locTemp(3) = jj;

                localizationResults(:,jj) = locTemp;

            end

            fprintf('Localization completed in %.2f seconds. \n',toc);
            fprintf('%.f percent of samples gave rise to a localization, for a total of %d localized points. \n',sum(~isnan(localizationResults(2,:))) / lengthData * 100,sum(~isnan(localizationResults(2,:))));


            %% Quality control?

            qualityControlThresh = 10; % mm;
            badInds = localizationResults(2,:) > qualityControlThresh;
            localizationResults(:,badInds) = [];

            %% Pack up

            self.sourceLocalizationResults.localizationResults = localizationResults;

            self.sourceLocalizationResults.paramStruct.marginError = marginError;
            self.sourceLocalizationResults.paramStruct.distanceThresh = distanceThresh;
            self.sourceLocalizationResults.paramStruct.perceivedActualCutoff = perceivedActualCutoff;

            self.sourceLocalizationResults.paramStruct.qualityControlThresh = qualityControlThresh;


        end

        function self = plotSurfFun(self)

            % p = inputParser;
            % addParameter(p,'timeWindow',0);
            % parse(p,varargin{:})
            % timeWindow = p.Results.timeWindow;


            %% Preamble

            currentEl = -90;

            %%%%
            % For seizure source info
            %%%%

            self = self.retrieveBraindata;
            myBp = self.braindata.myBp;
            myBd = self.braindata.myBd;

            isLeftInds = self.geodesic.isLeftInds;
            useLeft = logical(round(sum(isLeftInds) / length(isLeftInds)));

            %% Brain plot

            %%%%%
            % Basic setup
            %%%%%

            colormap(jet);
            figure
            myBd.ezplot(myBp,gca); % If we would not like to include the resection territory
            % plotResectionSurf(myStruct) % If we would like to include the resection territory
            if useLeft; view(-90,currentEl);
            else; view(90,currentEl);
            end
            myBp.camlights(5);
            ax = gca;

            %%  Grab and set up ROI results

            assert(~isempty(self.sourceLocalizationResults.localizationResults),'No data to plot.');
            self = self.locDataToRoi;
            roiResults = self.sourceLocalizationResults.roiResults;
            vMapAll = {roiResults.vertexMap};
            % clim = roiResults.clim;


            %% Find max color

            % The coloring is subject to the amount of ROI hits and so is difficult to
            % predict before mapping the data to ROIs. Therefore first we determine max
            % color, and then we plot.

            % This section is only needed for timeWindow ~= 0.

            maxColor = 0;
            for jj = 1:length(vMapAll)
                if isempty(vMapAll{jj}); continue; end

                vertexMap = vMapAll{jj};

                mapKeys = cell2mat(vertexMap.keys);
                VPerROI = cell(size(mapKeys));
                valPerROI = zeros(size(mapKeys));

                for ii = 1:length(mapKeys)
                    roiStruct = vertexMap(mapKeys(ii));
                    VPerROI{ii} = roiStruct.roi;
                    valPerROI(ii) = roiStruct.count;
                end
                isLeft = sign(mapKeys) == -1;

                if all(isLeft)
                    currentMax = findMaxColor(myBp, VPerROI, valPerROI,'surf','lh');
                elseif all(~isLeft)
                    currentMax = findMaxColor(myBp, VPerROI, valPerROI,'surf','rh');
                else
                    [isLeft,sortInds] = sort(isLeft,'descend');
                    VPerROI = VPerROI(sortInds);
                    valPerROI = valPerROI(sortInds);
                    rh_begin = find(~isLeft,1);

                    currentMax = findMaxColor(myBp, VPerROI, valPerROI,'rh_begin',rh_begin);

                end
                maxColor = max(maxColor,currentMax);
            end
            % fprintf('Max color is %.f. \n',maxColor);

            clim = [0 maxColor];

            %% Plot what we have

            lastPlot = false;

            for jj = 1:length(vMapAll)

                if ~isempty(vMapAll{jj})

                    vertexMap = vMapAll{jj};

                    mapKeys = cell2mat(vertexMap.keys);
                    VPerROI = cell(size(mapKeys));
                    valPerROI = zeros(size(mapKeys));

                    for ii = 1:length(mapKeys)
                        roiStruct = vertexMap(mapKeys(ii));
                        VPerROI{ii} = roiStruct.roi;
                        valPerROI(ii) = roiStruct.count;
                    end


                    isLeft = sign(mapKeys) == -1;

                    axes(ax);
                    if lastPlot; myBp.clearRegions(); end

                    if all(isLeft)
                        myBp.plotRegionsData(VPerROI, valPerROI,'surf','lh','clim',clim);
                    elseif all(~isLeft)
                        myBp.plotRegionsData(VPerROI, valPerROI,'surf','rh','clim',clim);
                    else
                        [isLeft,sortInds] = sort(isLeft,'descend');
                        VPerROI = VPerROI(sortInds);
                        valPerROI = valPerROI(sortInds);
                        rh_begin = find(~isLeft,1);

                        myBp.plotRegionsData(VPerROI, valPerROI,'rh_begin',rh_begin,'clim',clim);
                    end
                    fprintf('%d ROICs plotted; %d total ROI hits. \n',length(mapKeys),sum(valPerROI));
                    lastPlot = true;
                elseif lastPlot
                    axes(ax);
                    myBp.clearRegions();
                    lastPlot = false;
                end
                drawnow
            end


        end


        function self = locDataToRoi(self,varargin)

            p = inputParser;
            addParameter(p,'forceNew',false);
            % addParameter(p,'timeWindow',0);
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;
            % timeWindow = p.Results.timeWindow;

            % timeWindow appears to be deprecated (in the spikes version).
            % Maybe working in the seizure version?

            if ~forceNew && ~isempty(self.sourceLocalizationResults.roiResults); return; end

            self = self.localizationFunctionSpikes; % Needs to be done; no need to pass forceNew;

            self = self.retrieveBraindata;
            myBd = self.braindata.myBd;

            roiRadius = Inf;
            sumEmpty = 0;

            %%

            currentLocs = self.sourceLocalizationResults.localizationResults;
            currentLocs(:,isnan(currentLocs(1,:))) = [];

            assert(~isempty(currentLocs),'No localization.');


            currentVertices = currentLocs(1,:);

            uniqueVertex = unique(currentVertices);
            vertexMap = containers.Map('keyType','double','valueType','any');

            %%

            for vertIndex = uniqueVertex
                currentInds = find(currentVertices==vertIndex);
                numVertices = length(currentInds);

                isLeft = sign(vertIndex) == -1;

                if isLeft; surfString = 'lh'; currentSign = -1;
                else; surfString = 'rh'; currentSign = 1;
                end
                % Let's give left-sided vertices a negative sign.

                [roicUnique, ~, roiOutput] = myBd.vertex2ROI(abs(vertIndex),surfString,roiRadius);
                if isempty(roicUnique); sumEmpty = sumEmpty + numVertices; continue; end

                roicUnique = revertToVector(roicUnique);
                roicAll = roiOutput.ROIC_mesh_ndx;

                for roicIndex = roicUnique

                    %% Upload, second-wise

                    if isKey(vertexMap,roicIndex * currentSign)
                        roiStruct = vertexMap(roicIndex * currentSign);
                    else
                        roiStruct = struct;
                        roiStruct.count = 0;
                        roiStruct.roi = roiOutput.vertex(roicAll == roicIndex);
                        roiStruct.times = zeros(1,0);
                    end

                    roiStruct.count = roiStruct.count + numVertices;
                    roiStruct.times = [roiStruct.times currentLocs(3,currentInds)];

                    vertexMap(roicIndex * currentSign) = roiStruct;

                end
            end

            if sumEmpty; fprintf('%d unrecognized vertices found for %s.\n',sumEmpty,self.subj); end
            % fprintf('ROI data computed in %.1f seconds. \n',toc)

            roiResults = struct;
            roiResults.vertexMap = vertexMap;
            roiResults.timeWindow = 0;
            roiResults.roiRadius = roiRadius;

            self.sourceLocalizationResults.roiResults = roiResults;

        end



    end

    methods (Static = true)

        function [spkLeads,spkTimes] = removeDuplicates(spkLeads,spkTimes)


            seriesMap = containers.Map('keyType','char','valueType','any');
            for jj = 1:size(spkTimes,2)

                currentTimes = spkTimes(:,jj);
                currentTimes(isnan(currentTimes)) = [];
                for ii = 2:length(currentTimes)
                    currentTimes(ii) = sum(currentTimes(ii - 1:ii));
                end
                currentLeads = spkLeads(:,jj);

                % foundKey = false;
                for ii = 1:length(currentTimes)
                    currentKey = sprintf('%s_%d',currentLeads{ii},currentTimes(ii));

                    if isKey(seriesMap,currentKey)
                        seriesStruct = seriesMap(currentKey);
                        % foundKey = true;
                    else
                        seriesStruct = struct;
                        seriesStruct.indices = zeros(1,0);
                        seriesStruct.length = zeros(1,0);
                    end

                    seriesStruct.indices(end + 1) = jj;
                    seriesStruct.length(end + 1) = length(currentTimes);

                    seriesMap(currentKey) = seriesStruct;
                end

            end
            fullKeys = seriesMap.keys;

            removeIndices = zeros(1,0);
            for ii = 1:length(fullKeys)

                seriesStruct = seriesMap(fullKeys{ii});
                if isscalar(seriesStruct.indices); continue; end
                [~,ind] = max(seriesStruct.length);

                badIndices = seriesStruct.indices(1:length(seriesStruct.indices) ~= ind);
                removeIndices = [removeIndices badIndices];

            end
            removeIndices = unique(removeIndices);
            fprintf('%.2f%% of series removed as duplicates.\n',length(removeIndices) / size(spkLeads,2) * 100);
            spkLeads(:,removeIndices) = [];
            spkTimes(:,removeIndices) = [];
            % overallInds(removeIndices) = [];

            % warning('Remove duplicates done. Placing this warning so if we see it multiple times, we know this function is being called redundantly.');

        end

        function [rasters,waveforms] = removeVolCond_fromRaster(rasters,waveforms)

            listLocs = cell(1,size(rasters,2));

            if sum(rasters(:)) == 0; return; end


            volCondInds = find(sum(rasters,2) >= 2);
            if isempty(volCondInds); return; end
            fprintf('%.2f%% of single spike samples contained duplicate spikes and were removed on account of possible volume conduction.\n',full(length(volCondInds) / sum(any(rasters,2)) * 100));
            % markForDeletion = cell(1,size(rasters,2));

            for ii = 1:length(listLocs)
                listLocs{ii} = find(rasters(:,ii));
            end

            for ii = 1:length(listLocs)
                isBad = ismember(listLocs{ii},volCondInds);
                rasters(listLocs{ii}(isBad),ii) = false;
                waveforms{ii}(isBad,:) = [];
            end

        end



    end

end