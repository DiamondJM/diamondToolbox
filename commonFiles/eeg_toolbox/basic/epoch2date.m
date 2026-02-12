function dateOut = epoch2date(unix_epoch)
% EPOCH2DATE Translates (1000 times) a unix epoch to a MATLAB Serial date
%   Input:
%       unix_epoch - 1000 times Standard UNIX timestamp as number (microseconds since
%       01/01/1970).
%   Output:
%       dateOut - Serial date number (number of days since Jan 0, 0000).
%           This is the same type of date output by the DATENUM function.
%
%   Note: unix_epoch is not exactly a unix timestamp (which is in ms), but is expected to be 
%       1000 times a unix timestamp (ie microseconds). This is because this function is 
%       designed to read the timestamps of a session.log or eeg.log file, which are
%       1000*time.
%
%   See also DATENUM, DATESTR.
    if ~isnumeric(unix_epoch)
        error('epoch2date requires numeric input')
    end
    
    dt = datetime(unix_epoch / 1000, ...
                    'ConvertFrom', 'posixtime', ...
                    'TimeZone', 'local');
    
    dateOut = datenum(dt); 
end