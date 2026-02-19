function basePath = util_fill_params(pathInfoTable, basePath)
% It appears that this function works in the following way:
%
% paramKey  = {'PARAM_NAME', 'VAL'}
% inStr     = 'Has a [PARAM_NAME] parameter in brackets'
%
% outStr    = util_fill_params(paramKey, inStr)
%
% disp(outStr) >> 'Has a VAL parameter in brackets'
%
% This is basically string replacement

% convert input into proper table format
switch class(pathInfoTable)
    case 'table'
    case 'cell'
        pathInfoTable = cell2table(pathInfoTable,...
            'VariableNames',{'placeholder','val'});
    case 'struct' % this would also be a nice case to code up.
%                   % the only problem is that ... numerics make it hard.
        placeholder = fieldnames(pathInfoTable);
        val = struct2cell(pathInfoTable);
%         cellfun(@isstr,val);
        pathInfoTable = table(placeholder, val);
end
if isstr(basePath)
    % now... replace stuff
    for iRow = 1:height(pathInfoTable)
        fromStr = strcat('[', pathInfoTable{iRow,'placeholder'},']');
        toStr = pathInfoTable{iRow, 'val'};
        basePath = strrep(basePath, fromStr, toStr);
    end
    basePath = basePath{:};
else
end
% note: can technically do this all as cell array... may consider switching
% to tha instead, for less lines of code
return