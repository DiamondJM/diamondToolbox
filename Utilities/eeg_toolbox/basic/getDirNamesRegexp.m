function final_list = getDirNamesRegexp(dir_path, regular_exp, omit_list)
% Returns cell array of file and/or folder names you want based on a
% provided regular expression.
%
% final_list = getDirNamesRegexp(dir_path, regular_exp, omit_list)
% 
% Inputs:
%   dir_path (char): the path to the directory you want to return file and/or
%               folder names for
%   reg_exp (char): the regular expression of type char that will return the files and/or folders you
%               want
%   omit_list (cell array of chars): [optional] cell array of specific files and/or folders you definitely
%               do not want to include (probably want to keep the default)
%       default: omit_list = {'.','..','.DS_Store'};
%
% Outputs:
%   final_list (cell array of chars): cell array of the files and/or folders you want based on
%               the provided regular expression pattern
% 
% Examples:
%   getDirNamesRegexp(?/Volumes/Shares/FRNU/data/eeg?,?((NIH)|(TRE))\d{3}?) would give you all the subject folders in FRNU
%   getDirNamesRegexp(?/Volumes/Shares/FRNU/data/eeg/NIH066/raw?,?\d{6}_\d{4}?) would give you all session folders without the other files in the directory
% 
% 
% Created by Samantha Jackson
% 
% 10/2019 - Created by SJ
% 

if nargin < 3 % Set default omitting list
    omit_list = {'.','..','.DS_Store'};
end
if nargin < 2
    regular_exp = '.*';
end


% Start by returning the full directory
whole_dir = dir([dir_path filesep '*']); 

% Remove ".", "..", ".DS_Store"
mid_dir = whole_dir(~ismember({whole_dir.name},omit_list));

% Apply regular expression and return final_list
mid_list = {mid_dir.name};
final_list = mid_list(find(~cellfun(@isempty, regexp(mid_list,regular_exp))));

end

