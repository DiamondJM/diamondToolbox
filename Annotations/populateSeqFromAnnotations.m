function populateSeqFromAnnotations(sl, s, clipDetails, varargin)
% populateSeqFromAnnotations  Build sl.spikeDetectionResults.rasters from
% spreadsheetToSDTimes output.
%
%   populateSeqFromAnnotations(sl, s, clipDetails)
%   populateSeqFromAnnotations(..., 'subsensorLength', N)
%
%   sl          — sourceLocalizer object
%   s           — struct array from spreadsheetToSDTimes
%   clipDetails — struct from spreadsheetToSDTimes with fields:
%                   .clipStart    datetime
%                   .clipEnd      datetime
%                   .clipDuration duration
%
%   Builds a sparse (nSamples × nChans) raster matching the format used
%   by findSpikeTimes. Columns correspond to sl.chanNames; nSamples is
%   derived from clipDetails.clipDuration and sl.Fs.
%
%   Overwrites sl.spikeDetectionResults so only .rasters and .Fs remain
%   (any prior .waveforms / .paramStruct fields are dropped).

defaultLen = 3;
try
    defaultLen = sl.sourceLocalizationResults.paramStruct.subsensorLength;
catch
end

ip = inputParser;
ip.addParameter('subsensorLength', defaultLen, @isnumeric);
ip.parse(varargin{:});
subsensorLength = ip.Results.subsensorLength;

Fs        = sl.Fs;
chanNames = sl.chanNames;

assert(~isempty(chanNames), '[populateSeqFromAnnotations] sl.chanNames must be set.');
assert(isstruct(clipDetails) && isfield(clipDetails, 'clipDuration'), ...
    '[populateSeqFromAnnotations] clipDetails struct with clipDuration required.');

nSamples = round(seconds(clipDetails.clipDuration) * Fs);
nChans   = numel(chanNames);
raster   = sparse(nSamples, nChans);

% Filter to groups with at least subsensorLength electrodes
hasElec = arrayfun(@(g) numel(g.electrodes) >= subsensorLength, s);
s = s(hasElec);

nSpikes    = 0;
nUnmatched = 0;
nOutOfRange = 0;
for jj = 1:numel(s)
    elecs = s(jj).electrodes;

    % Sort electrodes by clip time
    clipSecs = arrayfun(@(e) seconds(e.clipTime), elecs);
    [clipSecs, ord] = sort(clipSecs, 'ascend');
    elecs = elecs(ord);

    % Convert to samples
    samps = round(clipSecs * Fs);

    for kk = 1:numel(elecs)
        idx = find(strcmp(chanNames, elecs(kk).name), 1, 'first');
        if isempty(idx)
            nUnmatched = nUnmatched + 1;
            continue
        end
        if samps(kk) < 1 || samps(kk) > nSamples
            nOutOfRange = nOutOfRange + 1;
            continue
        end
        raster(samps(kk), idx) = 1;
        nSpikes = nSpikes + 1;
    end
end

% Replace spikeDetectionResults — drop stale waveforms / paramStruct
sl.spikeDetectionResults = struct('rasters', raster, 'Fs', Fs);

fprintf(['[populateSeqFromAnnotations] %d spikes placed in %d × %d raster ' ...
    '(unmatched chan: %d, out-of-range: %d).\n'], ...
    nSpikes, nSamples, nChans, nUnmatched, nOutOfRange);
end
