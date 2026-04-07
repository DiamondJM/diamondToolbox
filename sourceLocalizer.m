classdef sourceLocalizer < handle

    properties

        subj

        chanNames
        Fs

        electrodeLocalizer

        braindata = struct('myBd',[],'myBp',[]);

        spikeDetectionResults = struct('rasters',[],'waveforms',[],'paramStruct',struct())
        seqResults = struct('seriesAll',[],'timesAll',[],'startEndTime',[]);
        sensors

        seizureProcessingResults = struct;

        sourceLocalizationResults = struct('localizationResults',[],'roiResults',[],'paramStruct',struct());
        localizationMode = 'spikes'; % or seizures

        deltaPosition

        timeSeries

    end

    properties (Hidden = true, Transient = true)

        subjFolder

    end

    properties (Transient = true)

        geodesic

    end

    properties (Hidden = true)
        rootFolder
    end

    methods

        %% Setup functions 

        function self = sourceLocalizer(subj, rootFolder, varargin)

            % Inputs:
            %   subj       - Subject name, e.g. 'NIH032'
            %   rootFolder - Path to root data folder.
            %       Expected folder structure:
            %           rootFolder/
            %               <subj>/
            %                   tal/
            %
            % Optional name-value:
            %   forceNewElectrodeLocalizer - If true, re-run electrode
            %       localization even if tal/leads.csv already exists.
            %       Default: false.
            %
            % Note:
            %   Channel names are prompted during the electrode localization
            %   create pipeline (just before the naming GUI). They can also
            %   be assigned directly (sl.chanNames) or loaded automatically
            %   from a time series file via sl.loadTimeSeries().
            %
            %   timeSeries and Fs are NOT constructor arguments. Either
            %   assign them directly after construction:
            %       sl.timeSeries = myData;   % [samples x channels]
            %       sl.Fs         = 1000;     % Hz
            %   or load from a file via dialog (.mat, .edf, .fif):
            %       sl.loadTimeSeries()
            %   Then call sl.localizationManager().

            p = inputParser;
            addParameter(p, 'forceNewElectrodeLocalizer', false);
            parse(p, varargin{:});
            forceNewElectrodeLocalizer = p.Results.forceNewElectrodeLocalizer;

            toolboxRoot = fileparts(mfilename('fullpath'));
            addpath(genpath(toolboxRoot));

            self.subj       = subj;
            self.rootFolder = rootFolder;
            self.chanNames  = {};

            %% Electrode localization

            self.electrodeLocalizer = electrodeLocalizer( ...
                subj, rootFolder, {}, 'forceNew', forceNewElectrodeLocalizer);

            % Propagate channel names resolved during electrode localization
            if ~isempty(self.electrodeLocalizer.chanNames)
                self.chanNames = self.electrodeLocalizer.chanNames;
            end

            %% Pull braindata
            self.retrieveBraindata;

        end

        %% Time series loading

        function loadTimeSeries(self)
            % Load a time series file via dialog and assign self.timeSeries,
            % self.Fs, and (when appropriate) self.chanNames.
            %
            % Supported formats:
            %   .mat — expects a 2-D numeric matrix [samples x channels] and
            %          a scalar Fs variable (searches for Fs/fs/srate/etc.).
            %          Falls back to an input dialog for Fs if not found.
            %   .edf — requires MATLAB R2023a+ Signal Processing Toolbox
            %          (edfread). Errors helpfully if unavailable.
            %   .fif — requires FieldTrip (ft_read_header / ft_read_data)
            %          on the MATLAB path. Errors helpfully if unavailable.
            %
            % Channel name resolution (chanNamesFromFile = names in the file):
            %   chanNames empty,  file has names   → set silently.
            %   chanNames set,    count matches     → leave existing.
            %   chanNames set,    count mismatch    → warn and overwrite
            %                                         (file header is
            %                                         authoritative for its
            %                                         own data).
            %   chanNames empty,  file has no names → error (ambiguous).
            %   chanNames set,    file has no names,
            %     count mismatch                    → error.
            %
            % Dismiss the dialog to abort without changes.

            [ts, Fs, chanNamesFromFile] = sourceLocalizer.loadTimeSeriesFromFile(self.chanNames);

            if isempty(ts)
                fprintf('[sourceLocalizer] No time series loaded.\n');
                return;
            end

            nCh = size(ts, 2);

            if ~isempty(chanNamesFromFile)
                if isempty(self.chanNames)
                    self.chanNames = chanNamesFromFile;
                elseif length(self.chanNames) ~= nCh
                    % Pre-set chanNames (e.g. from leads.csv) are the target set.
                    % Try to find each one in the file's channel list and subset ts.
                    [mask, idx] = ismember(strtrim(self.chanNames), strtrim(chanNamesFromFile));
                    if any(mask)
                        missing = self.chanNames(~mask);
                        if ~isempty(missing)
                            warning('[sourceLocalizer] %d channel(s) not found in EDF and will be missing: %s', ...
                                sum(~mask), strjoin(missing, ', '));
                        end
                        ts = ts(:, idx(mask));
                        self.chanNames = self.chanNames(mask);
                        nCh = size(ts, 2);
                        fprintf('[sourceLocalizer] Subsetted EDF to %d named channels.\n', nCh);
                    else
                        warning('[sourceLocalizer] chanNames (%d) does not match time series channels (%d) and no name overlap found. Overwriting with names from file.', ...
                            length(self.chanNames), nCh);
                        self.chanNames = chanNamesFromFile;
                    end
                end
                % else: existing chanNames match → leave alone
            else
                % File has no channel labels — prompt if we don't have them yet.
                if isempty(self.chanNames)
                    self.chanNames = sourceLocalizer.loadChanNamesFromFile();
                end
                if isempty(self.chanNames)
                    error('[sourceLocalizer] Channel names are required. Provide a channel names file or set sl.chanNames directly.');
                elseif length(self.chanNames) ~= nCh
                    error('[sourceLocalizer] chanNames has %d entries but time series has %d channels.', ...
                        length(self.chanNames), nCh);
                end
                % else: existing chanNames match → fine
            end

            self.timeSeries = ts;
            self.Fs         = Fs;
            fprintf('[sourceLocalizer] Loaded %d samples x %d channels at %.1f Hz.\n', ...
                size(ts,1), nCh, Fs);
        end

        %% Property setters

        function set.timeSeries(self, val)
            if ~isempty(self.chanNames) && ~isempty(val)
                assert(size(val, 2) == length(self.chanNames), ...
                    ['timeSeries must have %d columns (one per channel ' ...
                     'in chanNames), got %d.'], ...
                    length(self.chanNames), size(val, 2));
            end
            self.timeSeries = val;
        end

        function set.Fs(self, val)
            if ~isempty(val)
                assert(isnumeric(val) && isscalar(val) && val > 0, ...
                    'Fs must be a positive scalar.');
            end
            self.Fs = val;
        end

        function subjFolder = get.subjFolder(self)

            subjFolder = fullfile(self.rootFolder,self.subj);

        end


        function [myBd,myBp] = retrieveBraindata(self,varargin)

            %% Preamble

            p = inputParser;
            addParameter(p,'forceNew',false);
            addParameter(p,'saving',true);
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;
            saving = p.Results.saving;

            if ~forceNew && ~isempty(self.braindata.myBd)
                myBd = self.braindata.myBd; 
                myBp = self.braindata.myBp; 
                return; 
            end

            %%

            assert(exist('braindata2', 'file') == 2, 'braindata2 not found on path. Please check paths. Navigate to diamondToolbox and use addpath(genpath(pwd)).');
            assert(exist('brainplotter', 'file') == 2, 'brainplotter not found on path. Please check paths. Navigate to diamondToolbox and use addpath(genpath(pwd)).');

            subjDir = fullfile(self.subjFolder);

            fName = fullfile(subjDir,sprintf('brainData_%s.mat',self.subj));

            if exist(fName','file') == 2 && ~forceNew
                S = load(fName,'myBd','myBp');
                myBp = S.myBp;
                myBd = S.myBd;
            else

                fprintf('Creating braindata objects for %s.\n',self.subj);

                myBd = braindata2(self.subj,self.rootFolder);
                myBp = ez_get_plotter(myBd);

                if saving; save(fName,'myBd','myBp'); end

            end

            self.braindata.myBp = myBp;
            self.braindata.myBd = myBd;

        end

        %% Principal localization pipeline and helpers 

        function localizationManager(self,varargin)

            assert(~isempty(self.timeSeries), ...
                ['timeSeries must be set before running localizationManager. ' ...
                 'Assign it as: sl.timeSeries = yourData;']);
            assert(~isempty(self.Fs), ...
                ['Fs must be set before running localizationManager. ' ...
                 'Assign it as: sl.Fs = yourFs;']);

            p = inputParser;
            addParameter(p,'plotting',true);
            addParameter(p,'forceNew',false);
            parse(p,varargin{:})
            plotting = p.Results.plotting;
            forceNew = p.Results.forceNew;

            self.prepareDeltaPosition('forceNew',forceNew); % Differential handling for spikes and seizures

            %% Onto localization

            self.localizationFunction('forceNew',forceNew); % Spikes or seizure agnostic call
            self.locDataToRoi;

            %% Plot

            if ~plotting; return; end
            self.plotSurfFun;

            % self.plotDimensionsReducedWrapper; 

        end

        function prepareDeltaPosition(self,varargin)

            p = inputParser;
            addParameter(p,'forceNew',false);
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;

            if ~forceNew && ~isempty(self.deltaPosition); return; end

            switch self.localizationMode
                case 'spikes'
                    %% Preamble

                    % Define parameters
                    self.sourceLocalizationResults.paramStruct = struct(...
                        'propagationSpeed',300,... % mm / s
                        'sensorDistance',30,... mm
                        'subsensorLength',3,...
                        'distanceThresh', 30); % mm

                    % These are defaults and can be adjusted.
                    % propagationSpeed is in mm/s and refers to speed of wave
                    % propagation.

                    % sensorDistance is in mm. This is the distance sensors are apart
                    % from each other.

                    % subsensorLength = 3;
                    % How many electrodes at a time (at minimum) should participate in
                    % source localization.

                    % Here, since we're doing localization wtih spikes, 
                    % sensorDistance can basically be as far apart as we'd like, barring
                    % concerns for the distance at which the signal can travel.
                    % We do NOT have aliasing concerns like we did previously.
                    % Before, a large sensorDistance meant that we required very low
                    % frequencies, otherwise we'd have spatial aliasing.
                    % Here, interictal sourceLoc are rare (rather than periodic) events. We assume
                    % they arrive rarely, and so aliasing is not an issue. Since frequency of
                    % these events are taken to be ~0, there's no upper limit on inter-sensor
                    % distance (as dictated by frequency).

                    % Populate spikes
                    self.populateSpikes;

                    % Compute sequences
                    self.computeSequences;

                    self.getSensors;

                    %% Collect series

                    self.spikeSequenceToDeltaPosition;

                case 'seizure'

                    % Define parameters
                    self.sourceLocalizationResults.paramStruct = struct(...
                        'propagationSpeed',300,... % mm / s
                        'sensorDistance',12,... mm
                        'subsensorLength',4,...
                        'distanceThresh', 30); % mm

                    self.extractPhasePower;
                    self.postProcessPhase;

                    self.getSensors;

                    self.phaseToDeltaPosition;

            end

        end


        %% Ancillary functions, for IED localization 

        function populateSpikes(self,varargin)

            p = inputParser;
            addParameter(p,'forceNew',false);
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;

            if ~isempty(self.spikeDetectionResults.rasters) && ~forceNew; return; end

            assert(~isempty(self.timeSeries),'Time series empty.'); 
            %  Please assign timeSeries and Fs directly after construction:
            %       sl.timeSeries = myData;   % [samples x channels]
            %       sl.Fs         = 1000;     % Hz
            %   or load from a file via dialog (.mat, .edf, .fif):
            %       sl.loadTimeSeries()
            %   Then call sl.localizationManager().

            self.timeSeries = zscore(self.timeSeries);

            %% Define parameters

            self.spikeDetectionResults.paramStruct = struct(...
                'seqWin',0.1,... seconds
                'maxNegPeakWidth',0.05,... seconds
                'peakWin', 0.1, ... % seconds
                'ampScale',3,... z
                'trackPeaks',false,...
                'zThresh',3);

            fprintf('Calling spike detector with the following parameters. These are adjustable.\n');
            disp(self.spikeDetectionResults.paramStruct);

            self.findSpikeTimes;

        end

        function findSpikeTimes(self)

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


            fullRaster = sparse(tsDims(1),tsDims(2));
            waveformsMaster = cell(1,tsDims(2));

            parfor(kk=1:length(self.chanNames))

                warning('off','signal:findpeaks:largeMinPeakHeight');
                % Once for each worker 

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

        function computeSequences(self,varargin)

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

            startEndTime = nan(2,length(hasSeries));

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

                startEndTime(:,jj) = [totalWindows(1,hasSeries(jj));totalWindows(end,hasSeries(jj))];

            end

            seriesAll(:,all(isnan(timesAll))) = [];
            startEndTime(:,all(isnan(timesAll))) = [];
            timesAll(:,all(isnan(timesAll))) = [];
            % For sequences that became too short after accounting for duplicate
            % electrodes

            [seriesAll,timesAll,startEndTime] = self.removeDuplicates(seriesAll,timesAll,startEndTime);
            % For duplicate SEQUENCES, big difference

            self.seqResults.seriesAll = seriesAll;
            self.seqResults.timesAll = timesAll;
            self.seqResults.startEndTime = startEndTime;

        end


        function spikeSequenceToDeltaPosition(self)

            %% Preamble

            intervalMax = self.sourceLocalizationResults.paramStruct.sensorDistance / self.sourceLocalizationResults.paramStruct.propagationSpeed; % Seconds
            intervalMax = intervalMax * self.Fs;

            % subsensorLength = 3;

            % pairHits = zeros(length(chanNames));

            spkLeads = self.seqResults.seriesAll;
            spkTimes = self.seqResults.timesAll;

            self.getSensors;

            %% Prepare ingredients

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

            self.deltaPosition = seriesTable;

            % I think it's good to put this into the sensors struct,
            % because seriesTable should have the same height as
            % sensorInds.

        end

        %% Ancillary functions for seizure localization

        function extractPhasePower(self)

            assert(exist('seizureWindowFilt','file'),'This file should exist in diamondToolbox/Utilities/Filt. Please check paths.');

            thisFilter = @seizureWindowFilt;
            self.sourceLocalizationResults.paramStruct.filter = thisFilter;
            thisFilter = thisFilter();

            filtData = filtfilt(thisFilter.Numerator,1,self.timeSeries);
            hilbertData = hilbert(filtData);
            phase = angle(hilbertData);

            freqData = self.Fs * diff(unwrap(phase)) / (2 * pi);
            freqData = medfilt1(freqData,1000,'truncate');
            freqData(end + 1,:) = freqData(end,:);

            %% Process and pack

            pow = abs(hilbertData);
            realData = real(filtData);

            self.seizureProcessingResults.phase = phase;
            self.seizureProcessingResults.pow = pow;
            self.seizureProcessingResults.freq = freqData;
            self.seizureProcessingResults.real = realData;

        end

        function postProcessPhase(self)

            %%%%%%
            % Now post-process this data in such a way as to only pick the good data
            %%%%%

            % quantileVal = 0.10;

            phase = self.seizureProcessingResults.phase;
            pow = self.seizureProcessingResults.pow;

            retainChans = 14;

            for jj = 1:size(phase,1)

                currentPow = pow(jj,:);
                [~, ind] = sort(currentPow,'descend');

                lowPowInds = ind(retainChans + 1:end);

                phase(jj,lowPowInds) = nan;
            end

            pow(isnan(phase)) = nan;

            fprintf('%.f percent of all phase data will be used. \n',sum(~isnan(phase(:))) / numel(phase) * 100);

            self.seizureProcessingResults.phase = phase;
            self.seizureProcessingResults.pow = pow;

            self.sourceLocalizationResults.paramStruct.retainChans = retainChans;

        end

        function phaseToDeltaPosition(self)

            %% Preamble

            freqMaster = self.seizureProcessingResults.freq;
            phaseMaster = self.seizureProcessingResults.phase;
            sensorInds = self.sensors.sensorInds;

            freqMean = mean(freqMaster,2,'omitnan');

            % Let's find a bound on frequency difference informed by the doppler
            % effect.
            propagationSpeed = self.sourceLocalizationResults.paramStruct.propagationSpeed;
            expectedSourceSpeed = 20;
            acceptibleRatio = (propagationSpeed + expectedSourceSpeed) / (propagationSpeed - expectedSourceSpeed);

            lengthTs = size(self.seizureProcessingResults.phase,1);
            self.sourceLocalizationResults.paramStruct.acceptibleRatio = acceptibleRatio;

            propagationSpeed = self.sourceLocalizationResults.paramStruct.propagationSpeed;

            deltaPositionLocal = nan(size(sensorInds,1),lengthTs);

            for ii = 1:size(sensorInds,1)
                currentSensor = sensorInds(ii,:);

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
                % Center around 0 

                deltaPositionLocal(ii,validInds) = -(currentPhase(:,1) * propagationSpeed ./ (2 * pi * currentFreq(:,1)) ...
                    - currentPhase(:,2) * propagationSpeed ./ (2 * pi * currentFreq(:,2)));

                % Positive distance --> currentSensor(1) is CLOSER than sensor(2)

            end

            self.deltaPosition = deltaPositionLocal;

            fprintf('Finished computing deltaPosition, for seizure, subject %s. \n%.f percent of all computed phase differences used. \n',self.subj, sum(~isnan(deltaPositionLocal(:))) / numel(deltaPositionLocal) * 100);

        end

        %% Generic (spike or seizure) functions for setting up brain anatomy and electrode utilization

        function getSensors(self)

            self.loadGeodesic;
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

        function loadGeodesic(self,varargin)

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

            if ~isempty(self.geodesic) && ~isequal(length(self.chanNames),length(self.geodesic.leadNames)); forceNew = true; end
            if ~isempty(self.geodesic) && ~forceNew; return; end

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

            [myBd, myBp] = self.retrieveBraindata;

            fprintf('Building geodesic distances for %s. This can take up to 10 minutes. \n',self.subj);

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

        %% The actual localization function 
        % The workhorse 

        function localizationFunction(self,varargin)

            p = inputParser;
            addParameter(p,'forceNew',false);
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;

            %% Preamble

            if ~forceNew && ~isempty(self.sourceLocalizationResults.localizationResults); return; end

            %% Onto the function

            tic

            deltaPositionLocal = self.deltaPosition / self.Fs; % Seconds
            deltaPositionLocal = deltaPositionLocal * self.sourceLocalizationResults.paramStruct.propagationSpeed; % mm

            lengthData = size(deltaPositionLocal,2);
            sensorInds = self.sensors.sensorInds;

            myBp = self.braindata.myBp;

            self.loadGeodesic;

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

            % distanceThresh = 30;  % mm
            distanceThresh = self.sourceLocalizationResults.paramStruct.distanceThresh;
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
                currentSubsensor = find(abs(deltaPositionLocal(:,jj)) < distanceNorm * perceivedActualCutoff);
                % Perceived actual cutoff doesn't affect our source localization procedure itself. Rather, it
                % LIMITS the number of electrode pairs included in our procedure. The
                % reason is that, if deltaPosition(p,jj) is too close to distanceNorm(p),
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

                drSubsensor = deltaPositionLocal(currentSubsensor,jj);

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

            localizationResults(:,isnan(localizationResults(1,:))) = []; 

            qualityControlThresh = 10; % mm;
            badInds = localizationResults(2,:) > qualityControlThresh;
            localizationResults(:,badInds) = [];

            %% Pack up

            self.sourceLocalizationResults.localizationResults = localizationResults;

            self.sourceLocalizationResults.paramStruct.marginError = marginError;
            self.sourceLocalizationResults.paramStruct.distanceThresh = distanceThresh;
            self.sourceLocalizationResults.paramStruct.perceivedActualCutoff = perceivedActualCutoff;

            self.sourceLocalizationResults.paramStruct.qualityControlThresh = qualityControlThresh;

            self.sourceLocalizationResults.paramStruct.originalSize = lengthData; 

            self.sourceLocalizationResults.paramStruct.meta.subj = self.subj; 
            self.sourceLocalizationResults.paramStruct.meta.localizationMode = self.localizationMode;
            % Others? 

        end

        %% Plotting or post-localization analysis 

        function locDataToRoi(self,varargin)

            p = inputParser;
            addParameter(p,'forceNew',false);
            addParameter(p,'timeWindow',0); % Seconds
            parse(p,varargin{:})
            forceNew = p.Results.forceNew;
            timeWindow = p.Results.timeWindow;

            if ~forceNew && ~isempty(self.sourceLocalizationResults.roiResults); return; end

            if isequal(self.localizationMode,'spikes'); assert(~timeWindow);
            elseif isequal(self.localizationMode,'seizure') && ~timeWindow
                fprintf('\n%s\n', repmat('%', 1, 70));
                fprintf('NOTE: You are localizing with seizures, but timeWindow was passed\n');
                fprintf('as 0. This is perfectly reasonable, but if you''d like to call\n');
                fprintf('plotSurfFun and watch a video of the seizure as it unfolds over\n');
                fprintf('time, pass in timeWindow = 1 (in seconds) or any other value.\n');
                fprintf('%s\n\n', repmat('%', 1, 70));
            end

            self.localizationFunction; % Needs to be done; no need to pass forceNew;

            myBd = self.retrieveBraindata;

            roiRadius = Inf;
            sumEmpty = 0;

            thisLocalizationResults = self.sourceLocalizationResults.localizationResults;

            % lengthTs = size(thisLocalizationResults,2);
            lengthTs = self.sourceLocalizationResults.paramStruct.originalSize;

            if ~timeWindow
                timeWindowSamples = lengthTs; % Samples
            else; timeWindowSamples = timeWindow * self.Fs; % Samples
            end
            % Set timeWindow to lengthTs if 0, so that this just creates a
            % single step in timeBuffer.

            stepsBuffer = buffer(1:lengthTs,timeWindowSamples);
            vMapAll = cell(size(stepsBuffer,2),1);

            maxValAll = zeros(1,size(stepsBuffer,2));

            %%

            for jj = 1:size(stepsBuffer,2)

                maxVal = 0;

                currentTime = stepsBuffer(:,jj);
                currentTime = currentTime(logical(currentTime));

                currentLocs = thisLocalizationResults(:,ismember(thisLocalizationResults(3,:),currentTime));

                currentVertices = currentLocs(1,:);

                if isempty(currentVertices); continue; end
                uniqueVertex = unique(currentVertices);

                vertexMap = containers.Map('keyType','double','valueType','any');

                for vertIndex = uniqueVertex
                    currentInds = ismember(currentVertices,vertIndex);
                    numVertices = sum(currentInds);

                    isLeft = sign(vertIndex) == -1;

                    if isLeft; surfString = 'lh'; currentSign = -1;
                    else; surfString = 'rh'; currentSign = 1;
                    end
                    % Let's give left-sided vertices a negative sign.

                    [roicUnique, ~, roiOutput] = myBd.vertex2ROI(abs(vertIndex),surfString,roiRadius);

                    if isempty(roicUnique); sumEmpty = sumEmpty + numVertices; continue; end
                    roicUnique = roicUnique';
                    roiAll = roiOutput.ROIC_mesh_ndx;

                    for roicIndex = roicUnique

                        if isKey(vertexMap,roicIndex * currentSign)
                            roiStruct = vertexMap(roicIndex * currentSign);
                            roiStruct.count = roiStruct.count + numVertices;
                        else
                            roiStruct = struct;
                            roiStruct.count = numVertices;
                            roiStruct.roi = roiOutput.vertex(roiAll == roicIndex);
                        end

                        maxVal = max(maxVal,roiStruct.count);
                        vertexMap(roicIndex * currentSign) = roiStruct;
                    end

                end

                vMapAll{jj} = vertexMap;

                maxValAll(jj) = maxVal;

            end

            %%

            if sumEmpty; fprintf('%d unrecognized vertices found for %s.\n',sumEmpty,self.subj); end
            % fprintf('ROI data computed in %.1f seconds. \n',toc)

            roiResults = struct;
            roiResults.vertexMap = vMapAll;
            roiResults.timeWindow = timeWindow;
            roiResults.roiRadius = roiRadius;
            roiResults.maxValAll = maxValAll;

            self.sourceLocalizationResults.roiResults = roiResults;

        end

        function [maxRoic,roicTable] = findTopRoic(self)
            vertexMap = self.sourceLocalizationResults.roiResults.vertexMap{1};
            roicTable = zeros(length(vertexMap),3);
            if isempty(vertexMap); maxRoic = nan; return; end
            mapKeys = cell2mat(vertexMap.keys);
            for ii = 1:length(mapKeys)
                roiStruct = vertexMap(mapKeys(ii));
                roicTable(ii,1) = mapKeys(ii);
                roicTable(ii,2) = roiStruct.count;
                if isfield(roiStruct,'times'); roicTable(ii,3) = mean(roiStruct.times); end
            end
            [~,inds] = sort(roicTable(:,2),'descend');
            roicTable = roicTable(inds,:);
            maxRoic = roicTable(1,1);
        end

        function plotSurfFun(self)

            % p = inputParser;
            % addParameter(p,'timeWindow',0);
            % parse(p,varargin{:})
            % timeWindow = p.Results.timeWindow;

            %% Preamble

            currentEl = -90;

            %%%%
            % For seizure source info
            %%%%

            [myBd, myBp] = self.retrieveBraindata;

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
            self.locDataToRoi;
            roiResults = self.sourceLocalizationResults.roiResults;
            vertexMap = roiResults.vertexMap;
            % clim = roiResults.clim;

            maxColor = max(self.sourceLocalizationResults.roiResults.maxValAll); 
            colorLimits = [1 max(maxColor)];

            %% Plot what we have

            lastPlot = false;

            for jj = 1:length(vertexMap)

                if ~isempty(vertexMap{jj})

                    vertexMapLocal = vertexMap{jj};

                    mapKeys = cell2mat(vertexMapLocal.keys);
                    VPerROI = cell(size(mapKeys));
                    valPerROI = zeros(size(mapKeys));

                    for ii = 1:length(mapKeys)
                        roiStruct = vertexMapLocal(mapKeys(ii));
                        VPerROI{ii} = roiStruct.roi;
                        valPerROI(ii) = roiStruct.count;
                    end

                    isLeft = sign(mapKeys) == -1;

                    % axes(ax);
                    if lastPlot; myBp.clearRegions(); end

                    if all(isLeft)
                        myBp.plotRegionsData(VPerROI, valPerROI,'surf','lh','clim',colorLimits,'cmap',turbo);
                    elseif all(~isLeft)
                        myBp.plotRegionsData(VPerROI, valPerROI,'surf','rh','clim',colorLimits,'cmap',turbo);
                    else
                        [isLeft,sortInds] = sort(isLeft,'descend');
                        VPerROI = VPerROI(sortInds);
                        valPerROI = valPerROI(sortInds);
                        rh_begin = find(~isLeft,1);

                        myBp.plotRegionsData(VPerROI, valPerROI,'rh_begin',rh_begin,'clim',colorLimits,'cmap',turbo);
                    end
                    fprintf('%d ROICs plotted; %d total ROI hits. \n',length(mapKeys),sum(valPerROI));
                    lastPlot = true;
                elseif lastPlot
                    % axes(ax);
                    myBp.clearRegions();
                    lastPlot = false;
                end
                drawnow
            end

        end

        function plotDimensionsReducedWrapper(self,varargin)

            p = inputParser;
            addParameter(p,'colorMode','auto'); % timing; heatmap
            parse(p,varargin{:})
            colorMode = p.Results.colorMode;

            figure

            useMatrices = {[1 0 0; 0 1 0; 0 0 1],... axial 
                [1 0 0; 0 0 1; 0 1 0],... coronal 
                [0 0 1; 1 0 0; 0 1 0]}; %  sagittal

            if isequal(self.localizationMode,'seizure') && isequal(colorMode,'auto')
                fprintf('\n%s\n', repmat('%', 1, 70));
                fprintf('HINT: Try passing ''colorMode'',''heatmap'' to view results in\n');
                fprintf('terms of localization density, rather than in terms of timing.\n');
                fprintf('%s\n\n', repmat('%', 1, 70));
            end


            for ii = 1:length(useMatrices)

                subplot(1,length(useMatrices),ii) 
                self.plotDimensionReduced('coeff',useMatrices{ii},'colorMode',colorMode);

            end

        end

        function plotDimensionReduced(self,varargin)

            p = inputParser;
            addParameter(p,'colorMode','auto'); % timing; heatmap 
            addParameter(p,'coeff',eye(3)); % eye(3) for axial; [1 0 0; 0 0 1; 0 1 0] for coronal; [0 0 1; 1 0 0; 0 1 0] for sagittal 
            parse(p,varargin{:})
            colorMode = p.Results.colorMode;
            coeff = p.Results.coeff; 

            %% Preamble

            [~, myBp] = self.retrieveBraindata; 

            currentVertices = self.sourceLocalizationResults.localizationResults(1,:);
            locationsMaster = self.convertVerticesToLocations(myBp, currentVertices);

            assert(~isempty(locationsMaster),'No localization for %s.',self.subj); 

            %% Process dimension reduction.

            score = locationsMaster * coeff;
            % score = (locationsMaster - mu) * coeff;
            [a,b,c] = unique(currentVertices,'stable');
            score = score(b,:);

            % The above should index into a, where each element is the number of
            % appearances of a.

            % We can sort by size so we plot largest first

            % locsCounts = histcounts(categorical(c)); 
            % [~,sortCount] = sort(locsCounts,'ascend');
            % sortCount = randperm(size(score,1));
            % sortCount = 1:size(score,1);
            % sortCount = size(score,1):-1:1;
            % score = score(sortCount,:);
            % locsCounts = locsCounts(sortCount);
            % countVals = normalizeToBounds(locsCounts,[8 100]);  % Scaling

            inRange = squareform(pdist(locationsMaster));
            currentRad = 2;
            inRange = sum(inRange <= currentRad);
            inRange = inRange(b); 

            [inRange,sortCount] = sort(inRange,'descend'); 
            score = score(sortCount,:);

            % maxDensity = max(inRange) / (4/3 * pi * currentRad ^ 3);
            maxDensity = max(inRange);
            fprintf('Max density is %.2f.\n',maxDensity);

            % countVals = normalizeToBounds(inRange,[15 50],[1 108]); warning('Scaling.');
            countVals = normalizeToBounds(inRange,[15 50]);

            % locMode = self.sourceLocalizationResults.paramStruct.meta.localizationMode; 
            locMode = self.localizationMode; 

            if isequal(colorMode,'auto')

                switch locMode
                    case 'spikes'; colorMode = 'heatmap'; 
                    case 'seizure'; colorMode = 'timing'; 
                end
                if isequal(locMode,'seizure')
                    st = dbstack; 
                    if ~isequal(st(2).name,'sourceLocalizer.plotDimensionsReducedWrapper')
                        fprintf('\n%s\n', repmat('%', 1, 70));
                        fprintf('HINT: Try passing ''colorMode'',''heatmap'' to view results in\n');
                        fprintf('terms of localization density, rather than in terms of timing.\n');
                        fprintf('%s\n\n', repmat('%', 1, 70));
                    end
                end
            elseif isequal(colorMode,'timing')
                if isequal(locMode,'spikes'); warning('Spike localization results will be colored according to timing, but timing is probably meaningless, since you''re considering IED localization.'); end                % timing mode intended for seizures. 
            end

            switch colorMode

                case 'heatmap'
                    cmap = turbo;
                    colorVals = round(normalizeToBounds(inRange,[1 size(cmap,1)]));
                    colorVals = cmap(colorVals,:);

                case 'timing'
                    cmap = parula;

                    colorVals = zeros(size(a));

                    for jj = 1:length(a)
                        colorVals(jj) = mean(self.sourceLocalizationResults.localizationResults(3,self.sourceLocalizationResults.localizationResults(1,:)==a(jj)));
                    end
                    colorVals = round(normalizeToBounds(colorVals,[1 size(cmap,1)]));
                    colorVals = colorVals(sortCount);

                    colorVals = cmap(colorVals,:);

            end

            %% Plot

            % figure
            cla; hold on

            boundaryVert = [myBp.surfaces.pial_lh.vertices; myBp.surfaces.pial_rh.vertices];
            boundaryVert = boundaryVert * coeff;
            boundaryVert = boundaryVert(1:10:end,:);
            boundaryVert = double(boundaryVert(:,[1 2]));
            boundInds = boundary(boundaryVert(:,1),boundaryVert(:,2));
            plot(boundaryVert(boundInds,1),boundaryVert(boundInds,2),'k-','linewidth',2);
            
            for ii = 1:size(score,1)
                plot(score(ii,1),score(ii,2),'.', ...
                    'color',colorVals(ii,:),...
                    'MarkerSize',countVals(ii))
            end

            %% Formatting

            [d1, d0] = self.coeffToOrientation(coeff);

            axis equal
            set(gca,'DataAspectRatio',[1 1 1])
            % set(gca,'PlotBoxAspectRatio',[3 4 4])
            box on

            xticks(xlim); yticks(ylim);
            % xticklabels([d0(1,1) d1(1,1)]);
            % yticklabels([d0(2,1) d1(2,1)]);
            xticklabels({sprintf('%s',d0{1,1}), sprintf('%s',d1{1,1})});
            yticklabels({sprintf('%s',d0{2,1}), sprintf('%s',d1{2,1})});

        end


        % -----------------------------------------------------------------
        %% Staggered time series viewer
        % -----------------------------------------------------------------

        function plotTimeSeries(self, varargin)
            % Staggered multi-channel time series viewer with mini overview and drag scroll.
            %
            % Usage:
            %   sl.plotTimeSeries()
            %   sl.plotTimeSeries('winSec', 1)         % 1-sec window (default)
            %   sl.plotTimeSeries('stagger', 1)       % z-score unit spacing (default)
            %   sl.plotTimeSeries('showSeq', false)   % disable sequence overlay (default on)
            %
            % Navigation:
            %   Drag the blue window on the mini overview to scroll
            %   Left/Right arrow  — scroll 25% of window
            %   PageUp/PageDown   — scroll one full window
            %   Up/Down arrow     — scale amplitude up/down (stagger unchanged)
            %   Scale +/− buttons — same as Up/Down arrow
            %
            % Sequence overlay (when showSeq=true and seqResults is populated):
            %   Shaded band for each detected sequence (startEndTime window)
            %   Dots at the signal minimum for each participating channel

            ip = inputParser;
            ip.addParameter('winSec',  1,     @isnumeric);
            ip.addParameter('stagger', 1,     @isnumeric);
            ip.addParameter('showSeq', ~strcmp(self.localizationMode,'seizure'), @islogical);
            ip.parse(varargin{:});
            winSec  = ip.Results.winSec;
            stagger = ip.Results.stagger;
            showSeq = ip.Results.showSeq;

            assert(~isempty(self.timeSeries), '[plotTimeSeries] timeSeries is empty.');
            assert(~isempty(self.Fs),         '[plotTimeSeries] Fs is not set.');

            ts       = self.timeSeries;           % [samples x channels]
            nSamp    = size(ts, 1);
            nChan    = size(ts, 2);
            t        = (0:nSamp-1)' / self.Fs / 60;  % minutes
            totalSec = t(end);                        % now totalMin
            winSec   = min(winSec / 60, totalSec);    % convert input (seconds) → minutes

            % Z-score per channel
            mu = mean(ts, 1, 'omitnan');
            sd = std(ts,  0, 1, 'omitnan');
            sd(sd == 0) = 1;
            ts = (ts - mu) ./ sd;

            % Fixed stagger offsets (scale never changes these)
            offsets = ((nChan-1):-1:0) * stagger;   % 1 x nChan

            % Display decimation — cap main axes at 100,000 rendered points
            MAX_DISP_SAMP = 500000;
            dispDec = max(1, floor(nSamp / MAX_DISP_SAMP));
            if dispDec > 1
                fprintf('[plotTimeSeries] Display decimated %dx (%d → %d pts).\n', ...
                    dispDec, nSamp, floor(nSamp/dispDec));
            end
            t_disp   = t(1:dispDec:end);
            ts_base  = ts(1:dispDec:end, :);   % z-scored, decimated — no scale/stagger yet

            scale = 1;   % mutable amplitude scale; stagger is always offsets

            % Channel names
            if isempty(self.chanNames) || numel(self.chanNames) ~= nChan
                cNames = arrayfun(@(k) sprintf('Ch%d',k), 1:nChan, 'UniformOutput', false);
            else
                cNames = self.chanNames(:)';
            end
            [tickVals, tickIdx] = sort(offsets);
            tickLabels = cNames(tickIdx);

            % ── Sequence overlay data ────────────────────────────────────────
            % Pre-compute patch times and dot positions (scale-dependent y computed
            % at draw time and updated on changeScale).
            seqPatchTimes = zeros(0,2);   % [nSeq x 2] start/end in minutes
            dotTimes      = zeros(1,0);   % x positions (minutes)
            dotSigBase    = zeros(1,0);   % z-scored signal value at dot (no scale, no offset)
            dotOffsets    = zeros(1,0);   % channel stagger offset for each dot

            hasSeq = showSeq && isfield(self.seqResults,'startEndTime') && ...
                     ~isempty(self.seqResults.startEndTime);
            if hasSeq
                SET   = self.seqResults.startEndTime;   % 2 x nSeq  (samples)
                sAll  = self.seqResults.seriesAll;       % maxLen x nSeq
                tAll  = self.seqResults.timesAll;        % maxLen x nSeq  (samples, first abs + diffs)
                nSeq  = size(SET, 2);

                seqPatchTimes = SET' / self.Fs / 60;    % nSeq x 2, minutes

                for jj = 1:nSeq
                    col = sAll(:, jj);
                    col = col(~cellfun(@isempty, col));
                    nContacts = numel(col);
                    if nContacts == 0, continue; end

                    % Reconstruct absolute sample times (first entry is absolute,
                    % remaining are diffs)
                    absSamples = cumsum(tAll(1:nContacts, jj));

                    for kk = 1:nContacts
                        chanIdx = find(strcmp(cNames, col{kk}), 1);
                        if isempty(chanIdx), continue; end

                        tMin = absSamples(kk) / self.Fs / 60;
                        [~, dispIdx] = min(abs(t_disp - tMin));

                        dotTimes(end+1)   = t_disp(dispIdx);
                        dotSigBase(end+1) = ts_base(dispIdx, chanIdx);
                        dotOffsets(end+1) = offsets(chanIdx);
                    end
                end
            end

            % ── Figure ──────────────────────────────────────────────────────
            fig = figure('Name', sprintf('Time Series — %s', self.subj), ...
                'NumberTitle','off', 'Color','w', ...
                'WindowStyle','docked', ...
                'KeyPressFcn',           @onKeyPress, ...
                'WindowButtonDownFcn',   @onMouseDown, ...
                'WindowButtonMotionFcn', @onMouseMove, ...
                'WindowButtonUpFcn',     @onMouseUp, ...
                'WindowScrollWheelFcn',  @onScroll);

            % Main axes
            axMain = axes('Parent', fig, ...
                'Position', [0.09 0.26 0.89 0.71], ...
                'Color','w', 'Box','off', 'TickDir','out', 'FontSize',11);
            hold(axMain, 'on');

            % Sequence shaded windows (drawn first so traces appear on top)
            for jj = 1:size(seqPatchTimes, 1)
                t0 = seqPatchTimes(jj,1);  t1 = seqPatchTimes(jj,2);
                hp = patch(axMain, [t0 t1 t1 t0], [-1e4 -1e4 1e4 1e4], ...
                    [1 0.55 0.1], 'FaceAlpha',0.15, 'EdgeColor','none', 'HitTest','off');
                hp.YLimInclude = 'off';
            end

            hLines = plot(axMain, t_disp, ts_base * scale + offsets, ...
                'Color',[0.15 0.35 0.65], 'LineWidth',0.6);

            % Sequence dots
            hSeqDots = [];
            if ~isempty(dotTimes)
                hSeqDots = scatter(axMain, dotTimes, dotSigBase * scale + dotOffsets, ...
                    40, [0.85 0.2 0], 'filled', 'HitTest','off');
            end

            set(axMain, 'YTick',tickVals, 'YTickLabel',tickLabels);

            % Mini overview axes
            axMini = axes('Parent', fig, ...
                'Position', [0.09 0.07 0.89 0.10], ...
                'Color',[0.94 0.94 0.94], 'Box','off', ...
                'YColor','none', 'FontSize',9);
            hold(axMini, 'on');
            miniDec = max(1, floor(nSamp / 3000));
            plot(axMini, t(1:miniDec:end), ts(1:miniDec:end,:) + offsets, ...
                'Color',[0.55 0.55 0.55], 'LineWidth',0.3);
            xlim(axMini, [0 totalSec]);
            xlabel(axMini, 'Time (min)');

            % Window indicator patch
            drawnow limitrate;
            yl = ylim(axMini);
            hPatch = patch(axMini, ...
                [0 winSec winSec 0 0], [yl(1) yl(1) yl(2) yl(2) yl(1)], ...
                [0.20 0.45 0.85], 'FaceAlpha',0.25, ...
                'EdgeColor',[0.20 0.45 0.85], 'LineWidth',1.2, 'HitTest','off');

            % Scale buttons
            uicontrol('Parent',fig, 'Style','pushbutton', 'String','Scale +', ...
                'Units','normalized', 'Position',[0.09 0.01 0.07 0.04], ...
                'FontSize',11, 'Callback',@(~,~) changeScale(1.5));
            uicontrol('Parent',fig, 'Style','pushbutton', 'String','Scale −', ...
                'Units','normalized', 'Position',[0.17 0.01 0.07 0.04], ...
                'FontSize',11, 'Callback',@(~,~) changeScale(1/1.5));

            % Drag state
            isDragging   = false;
            dragSource   = '';    % 'mini' or 'main'
            dragStartX   = 0;
            dragStartWin = 0;

            setWindow(0);   % initialise XLim + XTick

            % ── Nested callbacks ────────────────────────────────────────────
            function setWindow(tStart)
                tStart = max(0, min(tStart, totalSec - winSec));
                xlim(axMain, [tStart, tStart + winSec]);
                yl2 = ylim(axMini);
                hPatch.XData = [tStart tStart+winSec tStart+winSec tStart tStart];
                hPatch.YData = [yl2(1) yl2(1) yl2(2) yl2(2) yl2(1)];
                % Ticks every 100 ms
                tickStep = 0.1 / 60;   % 100 ms in minutes
                ticks    = (ceil(tStart / tickStep) * tickStep) : tickStep : (tStart + winSec);
                set(axMain, 'XTick', ticks, ...
                    'XTickLabel', arrayfun(@(x) sprintf('%.3gs', x*60), ticks, 'UniformOutput', false));
            end

            function changeScale(factor)
                scale = scale * factor;
                ts_new = ts_base * scale + offsets;
                for k = 1:nChan
                    hLines(k).YData = ts_new(:, k);
                end
                if ~isempty(hSeqDots) && isvalid(hSeqDots)
                    hSeqDots.YData = dotSigBase * scale + dotOffsets;
                end
            end

            function onMouseDown(~,~)
                % Check mini axes first
                pt  = get(axMini, 'CurrentPoint');
                xl  = xlim(axMini);
                yl2 = ylim(axMini);
                if pt(1,1) >= xl(1) && pt(1,1) <= xl(2) && ...
                   pt(1,2) >= yl2(1) && pt(1,2) <= yl2(2)
                    setWindow(pt(1,1) - winSec/2);
                    isDragging   = true;
                    dragSource   = 'mini';
                    dragStartX   = pt(1,1);
                    dragStartWin = axMain.XLim(1);
                    return;
                end
                % Check main axes
                pt  = get(axMain, 'CurrentPoint');
                xl  = xlim(axMain);
                yl2 = ylim(axMain);
                if pt(1,1) >= xl(1) && pt(1,1) <= xl(2) && ...
                   pt(1,2) >= yl2(1) && pt(1,2) <= yl2(2)
                    isDragging   = true;
                    dragSource   = 'main';
                    dragStartX   = pt(1,1);
                    dragStartWin = axMain.XLim(1);
                end
            end

            function onMouseMove(~,~)
                if ~isDragging, return; end
                if strcmp(dragSource, 'mini')
                    pt = get(axMini, 'CurrentPoint');
                    setWindow(dragStartWin + (pt(1,1) - dragStartX));
                else
                    pt = get(axMain, 'CurrentPoint');
                    % Dragging right pulls content right → window moves left
                    setWindow(dragStartWin - (pt(1,1) - dragStartX));
                end
            end

            function onMouseUp(~,~)
                isDragging = false;
                dragSource = '';
            end

            function onScroll(~, evt)
                tNow = axMain.XLim(1);
                setWindow(tNow + winSec * 0.1 * evt.VerticalScrollCount);
            end

            function onKeyPress(~, evt)
                tNow = axMain.XLim(1);
                switch evt.Key
                    case 'rightarrow', setWindow(tNow + winSec * 0.25);
                    case 'leftarrow',  setWindow(tNow - winSec * 0.25);
                    case 'pagedown',   setWindow(tNow + winSec);
                    case 'pageup',     setWindow(tNow - winSec);
                    case 'uparrow',    changeScale(1.5);
                    case 'downarrow',  changeScale(1/1.5);
                end
            end

        end


    end

    methods (Static = true)

        function [spkLeads,spkTimes,startEndTime] = removeDuplicates(spkLeads,spkTimes,startEndTime)

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
            startEndTime(:,removeIndices) = [];
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

        function vert = convertLocationsToVertices(myBp,locations)

            vertL = myBp.surfaces.pial_lh.vertices;
            vertR = myBp.surfaces.pial_rh.vertices;

            stdVert = myBp.stdNumVert;

            vertAll = [vertL;vertR];

            [~,vert] = min(pdist2(vertAll,locations));

            isLeft = vert <= stdVert;
            vert(isLeft) = vert(isLeft) * -1;
            vert(~isLeft) = vert(~isLeft) - stdVert;

        end
        
        function locations = convertVerticesToLocations(myBp,vertices)

            % Negative numbers in vertices should indicate left hemisphere.
            % Vertices are bounded between 1 and 198812, the number of
            % vertices on a standard pial surface. 

            %% Preamble

            lVertices = myBp.surfaces.pial_lh.vertices;
            rVertices = myBp.surfaces.pial_rh.vertices;

            %% Let's distribute our data into a master map.

            locations = nan(length(vertices),3); 

            for ii = 1:length(vertices)
                isLeft = sign(vertices(ii)) == -1;
                if isLeft; currentVert = lVertices;
                else; currentVert = rVertices;
                end

                if isnan(vertices(ii)); continue; end

                locations(ii,:) = currentVert(abs(vertices(ii)),:);
            end

        end

        function [d1, d0] = coeffToOrientation(coeff)

            % Orientation is RAS
            oPos = {'Right','Anterior','Superior'};
            oNeg = {'Left','Posterior','Inferior'};
            d1 = cell(3,3);
            d0 = d1;

            if size(coeff,2) < 3
                coeffTemp = zeros(3);
                coeffTemp(:,1:size(coeff,2)) = coeff;
                coeff = coeffTemp;
            end

            orientationPC = eye(3) * coeff';

            for ii = 1:3
                [~,inds] = sort(abs(orientationPC(ii,:)),'descend');

                posInds = sign(orientationPC(ii,inds)) == 1;

                d1(ii, posInds) = oPos(inds(posInds));
                d1(ii, ~posInds) = oNeg(inds(~posInds));

                d0(ii, posInds) = oNeg(inds(posInds));
                d0(ii, ~posInds) = oPos(inds(~posInds));

            end

            % Each ROW of d provides the biggest, medium, and smallest influence of the
            % initial dimensions, into that row.

        end

    end

    methods (Static)

        function ensureFieldTrip()
            % Ensure FieldTrip is on the MATLAB path, auto-detecting the
            % installation if needed. Checks in order:
            %   1. ft_read_header already on path → already configured.
            %   2. ft_defaults on path → call it to finish setup.
            %   3. FIELDTRIP_HOME environment variable.
            %   4. Common install parent directories, scanning for
            %      subdirectories matching 'fieldtrip*' (handles versioned
            %      installs like fieldtrip-20260218); latest version wins.

            if exist('ft_read_header', 'file') == 2
                return;
            end

            if exist('ft_defaults', 'file') == 2
                addpath(fileparts(which('ft_defaults')));
                addpath(fullfile(fileparts(which('ft_defaults')), 'utilities'));
                ft_defaults;
                return;
            end

            home = char(java.lang.System.getProperty('user.home'));

            parentDirs = {
                getenv('FIELDTRIP_HOME');
                fullfile(home, 'Documents', 'MATLAB');
                fullfile(home, 'MATLAB');
                fullfile(home, 'fieldtrip');
                '/usr/local';
            };

            for i = 1:numel(parentDirs)
                parent = parentDirs{i};
                if isempty(parent) || exist(parent, 'dir') ~= 7
                    continue;
                end

                % Parent itself might be the FieldTrip root
                if exist(fullfile(parent, 'ft_defaults.m'), 'file') == 2
                    addpath(parent);
                    addpath(fullfile(parent, 'utilities'));
                    ft_defaults;
                    return;
                end

                % Scan for fieldtrip* subdirectories; sort descending so
                % latest version wins
                d = dir(fullfile(parent, 'fieldtrip*'));
                d = d([d.isdir]);
                if isempty(d), continue; end
                names = fliplr(sort({d.name}));
                for j = 1:numel(names)
                    candidate = fullfile(parent, names{j});
                    if exist(fullfile(candidate, 'ft_defaults.m'), 'file') == 2
                        addpath(candidate);
                        addpath(fullfile(candidate, 'utilities'));
                        ft_defaults;
                        return;
                    end
                end
            end

            error(['[sourceLocalizer] FieldTrip not found. Install from ' ...
                   'fieldtriptoolbox.org and either add it to your MATLAB ' ...
                   'path or set the FIELDTRIP_HOME environment variable.']);
        end

        function chanNames = loadChanNamesFromFile()
            % Prompt the user to select a file containing channel names.
            % Returns a column cell array of strings, or {} if dismissed.
            %
            % Supported formats:
            %   .mat — loads a variable containing a cell array of strings.
            %   .csv — looks for a column named name/chanName/label/channel/
            %          electrode (case-insensitive); falls back to first column
            %          for headerless files. Compatible with BIDS electrodes.tsv
            %          (rename .tsv to .csv or use the 'name' column directly).
            %   .fif — reads header only via FieldTrip ft_read_header;
            %          accepts both head-only files (e.g. *-head.fif) and
            %          full data files. Requires FieldTrip on the path.
            %   .edf — reads header only via edfinfo (R2023a+); no data loaded.

            chanNames = {};

            while true
                choice = dlgNonModal( ...
                    {'Channel names populate the electrode naming GUI so you can assign each detected contact to a named channel. You can also type names manually in the GUI if you prefer.', ...
                     '', ...
                     'Accepted formats:', ...
                     '  .mat  — MATLAB workspace containing a cell array of strings', ...
                     '  .csv  — name/chanName/label column, or first column if no header', ...
                     '  .fif  — MNE/FieldTrip file; channel names read from header', ...
                     '  .edf  — EDF/EDF+ file; channel names read from header only'}, ...
                    'Channel Names', ...
                    'Browse...', 'Skip');

                if isempty(choice) || strcmp(choice, 'Skip')
                    return;
                end

                [fname, fpath] = uigetfile( ...
                    {'*.mat;*.csv;*.fif;*.edf', 'Channel names file (*.mat, *.csv, *.fif, *.edf)'}, ...
                    'Select channel names file');

                if ~isequal(fname, 0), break; end
                % file picker cancelled — loop back to description dialog
            end

            fullPath = fullfile(fpath, fname);
            [~, ~, ext] = fileparts(fname);

            if strcmpi(ext, '.mat')
                S = load(fullPath);
                fields = fieldnames(S);
                if isscalar(fields)
                    chanNames = S.(fields{1});
                else
                    [idx, ok] = listdlg( ...
                        'ListString',   fields, ...
                        'SelectionMode','single', ...
                        'PromptString', 'Select variable containing channel names:');
                    if ok
                        chanNames = S.(fields{idx});
                    else
                        fprintf('No variable selected. chanNames will be empty.\n');
                        return;
                    end
                end

            elseif strcmpi(ext, '.csv')
                % Try to find a named column first (handles BIDS and similar
                % formats where the column is called 'name', 'label', etc.).
                T = readtable(fullPath, 'ReadVariableNames', true);
                knownCols = {'channame','channames','name','names', ...
                             'channel','channels','channel_name','channel_names', ...
                             'label','labels','electrode','electrodes'};
                colMatch = find(ismember(lower(T.Properties.VariableNames), knownCols), 1);
                if ~isempty(colMatch)
                    chanNames = T{:, colMatch};
                else
                    % No recognised header — treat as headerless, take col 1.
                    T = readtable(fullPath, 'ReadVariableNames', false);
                    chanNames = T{:, 1};
                end

            elseif strcmpi(ext, '.fif')
                sourceLocalizer.ensureFieldTrip();
                hdr       = ft_read_header(fullPath);
                chanNames = hdr.label;

            elseif strcmpi(ext, '.edf')
                info      = edfinfo(fullPath);
                chanNames = cellstr(info.SignalLabels);

            else
                error('Unsupported file type: %s. Use .mat, .csv, .fif, or .edf.', ext);
            end

            % Normalise to column cell array of char
            if isstring(chanNames)
                chanNames = cellstr(chanNames);
            elseif isnumeric(chanNames)
                chanNames = arrayfun(@num2str, chanNames, 'UniformOutput', false);
            end
            if ~iscell(chanNames)
                chanNames = cellstr(chanNames);
            end
            chanNames = chanNames(:);
        end

        function [timeSeries, Fs, chanNames] = loadTimeSeriesFromFile(selectedChannels)
            % Prompt the user to select a time series file and return its
            % contents. Returns empty arrays if the dialog is dismissed.
            %
            % selectedChannels — optional cell array of channel names to load.
            %   For EDF files, only those signals are read from disk (faster).
            %   Pass {} or omit to load all channels.
            %
            % Returns:
            %   timeSeries  [samples x channels] double
            %   Fs          scalar sampling frequency (Hz)
            %   chanNames   m x 1 cell array of channel name strings
            %               (empty cell if the file provides none)
            if nargin < 1, selectedChannels = {}; end

            timeSeries = [];
            Fs         = [];
            chanNames  = {};

            while true
                choice = dlgNonModal( ...
                    {'Select a time series file.', ...
                     '', ...
                     'Accepted formats:', ...
                     '  .mat  — MATLAB workspace with a [samples x channels] matrix', ...
                     '  .edf  — European Data Format', ...
                     '  .fif  — MNE/FieldTrip file'}, ...
                    'Load Time Series', 'Browse...', 'Cancel');
                if isempty(choice) || strcmp(choice, 'Cancel')
                    return;
                end
                [fname, fpath] = uigetfile( ...
                    {'*.mat;*.edf;*.fif', 'Time series files (*.mat, *.edf, *.fif)'; ...
                     '*.mat',             'MATLAB file (*.mat)'; ...
                     '*.edf',             'European Data Format (*.edf)'; ...
                     '*.fif',             'MNE FIF file (*.fif)'}, ...
                    'Select time series file');
                if ~isequal(fname, 0), break; end
                % cancelled file picker — loop back to description dialog
            end

            fullPath = fullfile(fpath, fname);
            [~, ~, ext] = fileparts(fname);

            switch lower(ext)
                case '.mat'
                    [timeSeries, Fs, chanNames] = sourceLocalizer.loadTsFromMat(fullPath);
                case '.edf'
                    [timeSeries, Fs, chanNames] = sourceLocalizer.loadTsFromEdf(fullPath, selectedChannels);
                case '.fif'
                    [timeSeries, Fs, chanNames] = sourceLocalizer.loadTsFromFif(fullPath);
                otherwise
                    error('[sourceLocalizer] Unsupported format: %s. Use .mat, .edf, or .fif.', ext);
            end
        end

        function [ts, Fs, chanNames] = loadTsFromMat(fullPath)
            % Load time series from a .mat file.
            % Searches for a 2-D numeric matrix and a scalar Fs variable.

            S      = load(fullPath);
            fnames = fieldnames(S);

            % Find 2-D numeric matrix candidates
            is2d  = cellfun(@(f) isnumeric(S.(f)) && ismatrix(S.(f)) && ~isscalar(S.(f)), fnames);
            candidates = fnames(is2d);

            if isempty(candidates)
                error('[sourceLocalizer] No 2-D numeric matrix found in %s.', fullPath);
            elseif isscalar(candidates)
                tsVar = candidates{1};
            else
                [idx, ok] = listdlg( ...
                    'ListString',   candidates, ...
                    'SelectionMode','single', ...
                    'PromptString', 'Select variable containing time series [samples x channels]:');
                if ~ok
                    ts = []; Fs = []; chanNames = {};
                    return;
                end
                tsVar = candidates{idx};
            end

            ts = double(S.(tsVar));

            % Locate Fs
            fsFields = {'Fs','fs','srate','SampleRate','sample_rate','samplingRate','Srate'};
            Fs = [];
            for i = 1:numel(fsFields)
                if isfield(S, fsFields{i}) && isscalar(S.(fsFields{i}))
                    Fs = double(S.(fsFields{i}));
                    break;
                end
            end
            if isempty(Fs)
                answer = inputdlg('Sampling frequency (Hz):', 'Enter Fs', 1, {'1000'});
                if isempty(answer)
                    error('[sourceLocalizer] Sampling frequency is required.');
                end
                Fs = str2double(answer{1});
            end

            % Locate channel names
            cnFields = {'chanNames','chan_names','chans','labels','channel_names','channels','label'};
            chanNames = {};
            for i = 1:numel(cnFields)
                if isfield(S, cnFields{i})
                    v = S.(cnFields{i});
                    if iscell(v) || isstring(v)
                        chanNames = cellstr(v(:));
                        break;
                    end
                end
            end
        end

        function [ts, Fs, chanNames] = loadTsFromEdf(fullPath, selectedChannels)
            % Load time series from an EDF file.
            % Requires MATLAB R2020b+ Signal Processing Toolbox (edfread/edfinfo).
            %
            % selectedChannels — optional cell array of signal label strings.
            %   When provided, only those channels are read from disk via
            %   edfread 'SelectedSignals', which is significantly faster for
            %   large EDF files with many unused channels.
            if nargin < 2, selectedChannels = {}; end

            assert(exist('edfread', 'file') == 2, ...
                ['edfread not found. EDF import requires MATLAB R2020b+ with the ' ...
                 'Signal Processing Toolbox.']);

            % Use edfinfo for metadata — returns a stable struct across versions
            info         = edfinfo(fullPath);
            Fs           = double(info.NumSamples(1)) / seconds(info.DataRecordDuration);
            allChanNames = cellstr(info.SignalLabels);

            % Resolve which channels to load
            if ~isempty(selectedChannels)
                [mask, ~] = ismember(strtrim(selectedChannels), strtrim(allChanNames));
                if ~all(mask)
                    missing = selectedChannels(~mask);
                    warning('[loadTsFromEdf] %d requested channel(s) not found in EDF: %s', ...
                        sum(~mask), strjoin(missing, ', '));
                end
                toLoad    = selectedChannels(mask);
                chanNames = toLoad(:);
                fprintf('[loadTsFromEdf] Loading %d of %d channels from EDF (subset via chanNames).\n', ...
                    numel(toLoad), numel(allChanNames));
            else
                toLoad    = {};
                chanNames = allChanNames;
            end

            % Use edfread for data — API varies by MATLAB version
            try
                % R2021b+: 'OutputFormat','array' returns [nSamples x nChan] double
                if ~isempty(toLoad)
                    data = edfread(fullPath, 'SelectedSignals', toLoad, 'OutputFormat', 'array');
                else
                    data = edfread(fullPath, 'OutputFormat', 'array');
                end
                ts = double(data);
            catch
                % R2020b fallback: returns timetable with one cell per record per channel
                if ~isempty(toLoad)
                    tbl = edfread(fullPath, 'SelectedSignals', toLoad);
                else
                    tbl = edfread(fullPath);
                end
                nChan = width(tbl);
                cols  = cell(1, nChan);
                for c = 1:nChan
                    cols{c} = vertcat(tbl{:, c}{:});
                end
                ts = double(horzcat(cols{:}));
            end
        end

        function [ts, Fs, chanNames] = loadTsFromFif(fullPath)
            % Load time series from an MNE FIF file.
            % Requires FieldTrip (ft_read_header / ft_read_data) on path.

            sourceLocalizer.ensureFieldTrip();

            hdr = ft_read_header(fullPath);
            dat = ft_read_data(fullPath, 'header', hdr);  % [channels x samples]

            ts        = dat';                % [samples x channels]
            Fs        = hdr.Fs;
            chanNames = hdr.label(:);        % column cell array
        end

    end

end