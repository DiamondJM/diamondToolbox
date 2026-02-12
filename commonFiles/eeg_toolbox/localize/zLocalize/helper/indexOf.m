function indices = indexOf(searchForMat, searchInVec)
% indexOf searches a vector for a set of values and returns their indices
%
% indices = indexOf(searchForMat, searchInVec)
%
% This is similar to find(A == x,1) except that this function can take a
% matrix for x.
%
% Works for strings too
% example:
%   classes = {'center', 'center_min', 'center_mid', 'center_max'};
%   class_vals = {'center' 'center_max' 'center_max' 'center_max' 'center_max'}
%   indexOf(class_vals, classes)
%   >  [1     4     4     4     4]     
%
% Inputs:
%   searchForMat - m x n matrix of values to search for
%   searchInVec - unique set of values to search
%
% Outputs:
%   indices - m x n matrix. indices(i,j) contains the index in which
%           searchForMat(i,j) was found in searchInVec, or 0 if not found
%
% See Also: find

 % TODO: we can probably speed this up by sorting the Mat first and slicing
 % the searchIn, then reverse sorting

    searchInVec = searchInVec(:);
    assert(length(unique(searchInVec)) == length(searchInVec), 'Search values must be unique');

    indices = zeros(size(searchForMat));
    
    siz = size(searchForMat);
    
    % experimental:
    
    
    parfor i = 1 : siz(1) * siz(2)
        val = searchForMat(i);
        indices(i) = indexOfScalar(searchInVec, val);
    end
end

function ndx = indexOfScalar(searchIn, val)
    % Adjust the equality for strings
    if iscellstr(val), val = char(val); end
    useChar = ischar(val);
    if useChar
        ndx = find(strcmp(searchIn, val), 1);
    else
        ndx = find(searchIn == val, 1);
    end
    
    if isempty(ndx), ndx = 0; end
end
