function s = addField(s, fieldname, vals)
% ADDFIELD adds a new field to a structure array with given values
%
% s = addField(s, fieldname, vals)
%
% INPUTS
%   s - structure to modify
%   fieldname - name of new field
%   vals (optional) - values to add to structure
%
% OUTPUTS
%   s - modified structure
%
%
%
% REVISION HISTORY
%   10/16 MST - Created
%   10/16 TCS - Showed MST everything


% default to empty cell if no vals
if nargin == 2
    vals = cell(size(s));
end

if isempty(s)
    s=struct();
end

% Numerics should be cells
if isnumeric(vals) || islogical(vals)
    c = num2cell(vals);
elseif iscell(vals)
    c = vals;
elseif numel(s) == 1
    s(:).(fieldname) = vals;
    return;
else
    error('Unrecognized type for given values. Expected numeric,logical, or cell');
end

if numel(s) == 1
    s = repmat(s,size(vals));
end
    

assert((length(c) == length(s)), 'Given values not of same length as given structure (%d vs %d)', length(c), length(s));



% Assignment

[s(:).(fieldname)] = deal(c{:});
    
end