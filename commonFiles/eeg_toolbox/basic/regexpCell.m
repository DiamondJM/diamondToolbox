function [cellMatch] = regexpCell(cellArray,expression)
%regexpCell : Super quick function just to avoid 2 lines of code every time
%you want to use regexp on cell arrays
% 
% Inputs:
%       cellArray   = your cell array of strings
%       expression  = your regular expression as you woul duse with regexp
%
% Output:
%       cellMatch   = your un-nested cell array containing the matches deom
%                     expression
%
% Created by Samantha N Jackson (7/9/2020)
%

nestedCell = regexp(cellArray,expression,'match');

cellMatch = [nestedCell{:}];

end

