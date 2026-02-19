function [newF, info] = newerFile(file1,file2)
% 
% Compare the modification date of two files and return the name of the newer one.
%
% Modification of functions found in eegPrepAndAlign: isOldJacksheet and getFileModDate
% Returns the newer file as newF. If they were modified at the same time, newF is empty ''
% 
% inputs:
%       file1: string to the path of the file
%               ex) file1 = '/Volumes/SeagateBackupPlusDrive/local56/56PROC/eeg_processing/.update/eeg/NIH062/docs/element_info.csv'
%       file2: string to the path of the file
%
% outputs:
%       newF: A string specifying the file path of either file1 or file2, whichever is newer, OR an empty
%               string if they were modified at the same time
%       info: A cell array of size 2x1 containing the modification date and file path of file1 and file2
%
%
% Created by Samantha N. Jackson 2/3/2020
%

% file1:
[~,result1] = unix(sprintf('stat -f %%m "%s"', file1)); %maybe switch to %%c ? ; SJ 9/19/19 - change %s to "%s"
result1 = strsplit(result1);
date1 = epoch2date(1000 * str2double(result1{1}));

% file2:
[~,result2] = unix(sprintf('stat -f %%m "%s"', file2)); %maybe switch to %%c ? ; SJ 9/19/19 - change %s to "%s"
result2 = strsplit(result2);
date2 = epoch2date(1000 * str2double(result2{1}));

if date2 > date1 %file2 is newer
    newF = file2;    
elseif date2 == date1
    newF = '';
else %file1 is newer
    newF = file1;
end

if nargout > 1
    info = {[datestr(date1) ': ' file1], [datestr(date2) ': ' file2]}';

end

end