function [grouped, group_ndxs] = groupBySubstring(longerStrings, groupBy)
% GROUPBYSUBSTRING groups a string array into multiple string arrays by substring
%
% [grouped, group_ndxs] = groupBySubstring(longerStrings, groupBy)
%
% This is not generic! It only works for tag names/ channel names
% It depends upon a channel name being a tag name plus a number like
% TG and TG2
%
%
% Inputs:
%   longerStrings: cell array of strings to group
%   groupBy: cell array of shorter strings that determine group
%
% Output:
%   grouped: cell array, each element is cell array of longerString members
%   group_ndxs: cell array, i'th element is an array of indexes for the i'th group
%
% Notes: case-insentitive
    
    groupBy = cellstr(groupBy);
    n = length(groupBy);
    grouped = cell(1, n);
    group_ndxs = cell(1, n);
    
    for i = 1 : n
        substring = groupBy{i};
        matchMask = strcmp(util_split_stringnum(longerStrings), substring);
        group = longerStrings(matchMask);
        grouped(i) = {group};
        group_ndxs{i} = find(matchMask);
    end
    

end