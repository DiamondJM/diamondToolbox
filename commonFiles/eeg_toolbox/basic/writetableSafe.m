function writetableSafe(t, filename, varargin)
% WRITETABLESAFE calls writetable and catches common errors
%
% Error-handling: 
%   comma's within a cell - replaces comma's with spaces so that MATLAB_2015b can readtable without issue
%   
    c = table2cell(t);
    for i = 1:numel(c)
        if ischar(c{i})
            c{i} = strrep(c{i}, ',', ' ');
        end
    end
    t = cell2table(c, 'VariableNames', t.Properties.VariableNames);
    writetable(t, filename, varargin{:});
end