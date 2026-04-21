function dt = edfStartTime(edfPath)
% edfStartTime  Read recording start datetime from an EDF file header.
%
%   dt = edfStartTime(edfPath)
%
%   Reads the fixed-width EDF header directly (no toolbox required).
%   EDF standard: date bytes 168-175 'dd.mm.yy', time bytes 176-183 'hh.mm.ss'.

fid = fopen(edfPath, 'r', 'n', 'US-ASCII');
assert(fid ~= -1, '[edfStartTime] Cannot open file: %s', edfPath);
fseek(fid, 168, 'bof');
dateStr = fread(fid, 8, '*char')';
timeStr = fread(fid, 8, '*char')';
fclose(fid);

% Parse 'dd.mm.yy' and 'hh.mm.ss'
try
    dt = datetime(sprintf('%s %s', dateStr, timeStr), ...
        'InputFormat', 'dd.MM.yy HH.mm.ss');
catch
    error('[edfStartTime] Could not parse EDF date/time: "%s" "%s"', dateStr, timeStr);
end
end
