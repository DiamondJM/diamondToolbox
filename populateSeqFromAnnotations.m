function populateSeqFromAnnotations(sl, s, varargin)
% populateSeqFromAnnotations  Fill sl.seqResults from spreadsheetToSDTimes output.
%
%   populateSeqFromAnnotations(sl, s)
%   populateSeqFromAnnotations(sl, s, 'subsensorLength', N)
%
%   sl  — sourceLocalizer object
%   s   — struct array from spreadsheetToSDTimes
%
%   Populates sl.seqResults.seriesAll, timesAll, and startEndTime.
%   Only retains SD groups whose electrode count equals subsensorLength.
%
%   timesAll format (matches computeSequences convention):
%     row 1     : absolute sample of first electrode detection from clip onset
%     rows 2..N : sample offsets from row 1 (cumsum reconstructs absolute positions)

% Resolve default subsensorLength from paramStruct, fallback to 3
defaultLen = 3;
try
    defaultLen = sl.sourceLocalizationResults.paramStruct.subsensorLength;
catch
end

ip = inputParser;
ip.addParameter('subsensorLength', defaultLen, @isnumeric);
ip.parse(varargin{:});
subsensorLength = ip.Results.subsensorLength;

Fs = sl.Fs;

% Filter to groups with exactly subsensorLength electrodes
hasElec = arrayfun(@(g) numel(g.electrodes) == subsensorLength, s);
s = s(hasElec);
nSeq = numel(s);

if nSeq == 0
    warning('[populateSeqFromAnnotations] No SD groups with electrode entries found.');
    return
end

% Determine max electrodes per group for pre-allocation
maxLen = max(arrayfun(@(g) numel(g.electrodes), s));

seriesAll   = cell(maxLen, nSeq);
seriesAll(:) = {''};
timesAll    = nan(maxLen, nSeq);
startEndTime = nan(2, nSeq);

for jj = 1:nSeq
    elecs = s(jj).electrodes;
    nElec = numel(elecs);

    % Sort electrodes by clip time
    clipSecs = arrayfun(@(e) seconds(e.clipTime), elecs);
    [clipSecs, ord] = sort(clipSecs, 'ascend');
    elecs = elecs(ord);

    % Convert to samples
    samps = round(clipSecs * Fs);

    % seriesAll: electrode names in time order
    for kk = 1:nElec
        seriesAll{kk, jj} = elecs(kk).name;
    end

    % timesAll: first row absolute, rest are diffs
    timesAll(1, jj) = samps(1);
    if nElec > 1
        timesAll(2:nElec, jj) = diff(samps);
    end

    % startEndTime: group onset and onset+duration in samples
    groupSamp = round(seconds(s(jj).groupClipTime) * Fs);
    startEndTime(1, jj) = groupSamp;
    startEndTime(2, jj) = groupSamp + round(seconds(s(jj).durationOverall) * Fs);
end

sl.seqResults.seriesAll   = seriesAll;
sl.seqResults.timesAll    = timesAll;
sl.seqResults.startEndTime = startEndTime;

fprintf('[populateSeqFromAnnotations] %d sequences loaded.\n', nSeq);
end
