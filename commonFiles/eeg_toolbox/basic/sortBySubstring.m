function channelOrder = sortBySubstring(channels, tags)
% SORTBYSUBSTRING sorts an cell array of strings by the order of left substrings given in a smaller array.
%
% INPUT
%  channels - Strings to sort
%  tags     - shorter strings, "key" to sort by
%
% OUTPUT
%  channelOrder contains the order of sorted channels
%
% NOTES
%  - Sorts on longest match to a tag name
%  - Within substring matches, sorts numerically
%
%
% Example:
% channels = {'BAG1','BAG2','TAG3','TAG20','TAG1'}
% tags = {'TAG','BAG'}
% order = sortBySubstring(channels,tags)
% disp(channels(order)) --> {'TAG1','TAG3','TAG20','BAG1','BAG2'}
%
%
% Used to sort jacksheetMaster by order given in element_info
%
% REVISION HISTORY
%   08/16 MST Created
%   10/16 MST Added numerical sort within tag groups

channelOrder = zeros(length(channels), 1);
channelNdxByTag = zeros(length(tags), length(channels));
numTagMatches = zeros(length(tags));
assert(length(tags)==length(unique(tags)), 'Tags must be unique');

% if ~exist('presortChannels', 'var') || presortChannels
%    channels = sort(channels);
% end

% loop chanels in order
for channelNdx = 1 : length(channels)
    channel = channels{channelNdx};
    hiMatchLen = 0;
    matchTagNdx = 0;

    % loop tags in order
    for j = 1 : length(tags)
        tag = tags{j};
        tagLen = length(tag);

        % find the longest string match for this channel amongst tags
        if tagLen > hiMatchLen && strncmp(channel, tag, tagLen)
            hiMatchLen = tagLen;
            matchTagNdx = j;
        end
    end

    if matchTagNdx == 0, continue; end;

    % for each tag, record in order the channels that matched (by ndx)
    numTagMatches(matchTagNdx) = numTagMatches(matchTagNdx) + 1;
    channelNdxByTag(matchTagNdx, numTagMatches(matchTagNdx)) = channelNdx;

end

% The order we want is left-to-right by rows, so take nonzero's (gets
% columns) of transpose
channelOrder = nonzeros(channelNdxByTag');


% Now split and order by numbers within each tag
channelsSorted = channels(channelOrder);
[tags_sorted,nums] = util_split_stringnum(channelsSorted);

for i = 1 : length(tags)
    match_ndx = strcmpi(tags_sorted, tags{i});
    [~,subOrder] = sort(nums(match_ndx));
    tempOrder = channelOrder(match_ndx);
    channelOrder(match_ndx) = tempOrder(subOrder);
end


end % end function