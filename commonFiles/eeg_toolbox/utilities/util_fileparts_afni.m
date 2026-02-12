function [filedir, filename] = util_fileparts_afni(filepath, varargin)
% Desc:
%   for the afni file path of 
%   [filedir]/[filename]+orig
% Inputs:
%   filepath
% Outputs:
%   filedir
%   filename without '+orig'


[filedir, filename, ext] = fileparts(filepath);
filename = [filename ext];

if ~isempty(varargin)
    if strcmp(varargin{1},'full')
        % do nothing
    end
else
    filename = strrep(filename, '+orig', '');
end
return
