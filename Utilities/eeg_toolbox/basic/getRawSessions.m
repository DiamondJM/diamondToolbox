function [rawSessions, isStim] = getRawSessions(rawPath)
% getRawSessions returns all of the session folder names in eeg/NIH000/raw, including those in the
% STIM_MAP folder.
% Note: The sessions are sorted in chronological order
%
% Inputs:
%       rawPath           = (string) Path to a subject's raw directory where the sessions are
%                           ex) rootEEGdir = '/Volumes/SeagateBackupPlusDrive/data/eeg/NIH076/raw'
%
% Outputs:
%       rawSessions       = (cell) Array of all of the session folder names in 'raw' and 'raw/STIM_MAP'
%       isStim (Optional) = (logical) Array corresponding to the indices in rawSessions of true for session 
%                           from STIM_MAP and false for directly in raw folder
%
% Created by Samantha Jackson, 3/10/2020
% 
%
%
% 3/10/2020 - Created by SJ
% 


if ~exist(rawPath,'dir')
    fprintf('%s\n',['ERROR!!!! ' rawPath ' does not exist! Ask SJ']);
    keyboard
else
    % Get all of the session folders directly in 'raw':
    rawSessions = getDirNamesRegexp(rawPath,'^\d{6}_\d{4}.*');
    isStim = zeros(1,numel(rawSessions));
    % Now get any in STIM_MAP
    stimpath = fullfile(rawPath,'STIM_MAP');
    if ~exist(stimpath,'dir')
        fprintf('%s\n','ERROR!!! No STIM_MAP folder in raw directory? Why?');
        keyboard
    else
        stimSessions = getDirNamesRegexp(stimpath,'^\d{6}_\d{4}.*');
        rawSessions = [rawSessions, stimSessions];
        isStim = [isStim ones(1,numel(stimSessions))];
    end
    % Now sort them in chronological order
    [rawSessions, sortIdx] = sort(rawSessions);
    isStim = logical(isStim(sortIdx));
end
end

