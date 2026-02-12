function [tagName, absElmtIdx] = util_split_stringnum(tagIdx)
    % Desc: splits combined string and number into the string and number
    % Inputs:
    %   array of strings, combined string and number
    % Outputs:
    %   string
    %   number
    % Exampls:
    %   [str num] = util_split_stringnum('G15')
    %   str = 'G'
    %   num = 15
    % add optional funcionality to work on tables as well, iterating through
    % rows

    tagIdx = cellstr(tagIdx);

    tagName = cell(length(tagIdx), 1);
    absElmtIdx = zeros(length(tagIdx), 1);
    checkchar = @(myStr) ~isempty(myStr) && all(ismember(myStr,'0123456789'));

    for iTag = 1:length(tagIdx)
        myTagIdx = tagIdx{iTag};
        whichChar = arrayfun(@(x) checkchar(x), myTagIdx);
        ndx = myTagIdx(whichChar);
        if isempty(ndx)
            ndx = NaN; 
        else
            ndx = str2num(ndx);
        end
        tagName{iTag} = myTagIdx(~whichChar);
        absElmtIdx(iTag) = ndx;
    end

end