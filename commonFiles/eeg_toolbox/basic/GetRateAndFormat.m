function [samplerate,nBytes,dataformat,gain,sourceType] = GetRateAndFormat(event)
%GETRATEANDFORMAT - Get the samplerate, gain, and format of eeg data.
%
% function [samplerate,nBytes,dataformat,gain] = GetRateAndFormat(event)
%
% Input:
%   event - a struct that contains a field named 'eegfile' which
%           points to a folder that contains a params.txt file.
%               OR
%           a string which is the path to a folder that contains a
%           params.txt file
%
% Revision History
%   ?? - Unkown
%   09/17 MST - Works with non-session specific path (e.g. */eeg.processed/)
%   01/19 JW  - add sourceType.. string sayign what recording device was used

if ischar(event)
    % event is actually a path to an EEG file
    path = event;
else
    path = event.eegfile;
end


sessionParamsPath = path; %- possibly sending in full path to params.txt


%- want it to be a complete path to params.txt.
if ~contains(sessionParamsPath,'params.txt'),
    if exist(sessionParamsPath,'file'),
        [sessionParamsPath,NAME,EXT] = fileparts(sessionParamsPath);
    end
    sessionParamsPath = fullfile(sessionParamsPath,'params.txt');
end


% *MST 09/17 - code necessary to work-around original trimming from above block
% Tries three paths:
%   1) If session string given, try the general [session]/../eeg.*/params.txt
%   2) Try the actual path the user gave us: [path]/params.txt
%   3) If session string given, try the noreref general: [session]/../eeg.noreref/params.txt
file = fopen(sessionParamsPath,'rt'); 
if( file == -1 )
    % *MST 09/17 - try the actual path the user gave us
    badParamsPath = sessionParamsPath;
    sessionParamsPath = fullfile(path, 'params.txt');
    file = fopen(sessionParamsPath,'rt'); 
    if file == -1
        % *MST 09/17 - try the noreref folder
        badParamsPath2 = sessionParamsPath;
        sessionParamsPath = fullfile(path, '../..', 'eeg.noreref', 'params.txt');
        file = fopen(sessionParamsPath,'rt'); 
        if file == -1
            error('params not found. tried %s, %s, and %s', badParamsPath, badParamsPath2, sessionParamsPath);
        end
    end
end
pathUsed = fileparts(sessionParamsPath);
params = eegparams({'samplerate','gain','dataformat'},file);


samplerate = params{1};
if( isempty(samplerate) )
    % EH - can't do anything if the sample rate isn't present (no default)
    gain = [];
    dataformat = '';
    nBytes = [];
    fprintf('\nERROR: params file not found, or doesnt contain samplerate!\n');
    keyboard
    return;
end

if( ~isempty(params{2}) )
    gain = params{2};
else
    gain = 1.0;
end

if( ~isempty(params{3}) )
    dataformat = params{3};
else
    dataformat = 'short';
end

switch dataformat
    case {'short','int16'}
        nBytes = 2;
    case {'single','int32'}
        nBytes = 4;
    case {'double','int64'}
        nBytes = 8;
    otherwise
        error('BAD DATA FORMAT!');
end

%- also open the sourceType.txt file found in the same place
sourceType = 'nihon khoden raw';
sourceFile = fullfile(pathUsed,'sourcetype.txt');
if exist(sourceFile,'file'),
    fid = fopen(sourceFile,'r');
    sourceType = fgetl(fid);
    fclose(fid);
end

