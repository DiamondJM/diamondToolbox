function [anovaTable,fullTable] = rAnovaWrapper(resultsAll,varargin)

p = inputParser;
addParameter(p,'varNames',strsplit(string(num2str(1:size(resultsAll,2)))));
parse(p,varargin{:})
varNames = p.Results.varNames;

rTable = array2table(resultsAll); 

withinDesign = table(varNames','VariableNames',{'State'});
withinDesign.State = categorical(withinDesign.State);

conditionString = sprintf('%s-%s ~ 1',rTable.Properties.VariableNames{1},rTable.Properties.VariableNames{end}); 

rm = fitrm(rTable,conditionString,'WithinDesign',withinDesign);
fullTable = ranova(rm,'WithinModel','State');
anovaTable = simpleAnovaTable(fullTable,'Measure (units)');
disp(anovaTable)

% c = multcompare(rm,'State','ComparisonType','bonferroni');
% disp(c); 
% Skipping multcompare because MATLAB doesn't seem to offer 1)
% Holm-Bonferroni correction, or 2) pooling error term for repeated
% measures factors. For these reasons, I used JASP for multiple
% comparisons.

end

function [s] = simpleAnovaTable(AT, dvName)
c = table2cell(AT);
% remove erroneous entries in F and p columns 
for i=1:size(c,1)       
        if c{i,4} == 1
            c(i,4) = {''};
        end
        if c{i,5} == .5
            c(i,5) = {''};
        end
end
% use conventional labels in Effect column
effect = AT.Properties.RowNames;
for i=1:length(effect)
    tmp = effect{i};
    tmp = erase(tmp, '(Intercept):');
    tmp = strrep(tmp, 'Error', 'Participant');
    effect(i) = {tmp}; 
end
% determine the required width of the table
fieldWidth1 = max(cellfun('length', effect)); % width of Effect column
fieldWidth2 = 57; % field needed for df, SS, MS, F, and p columns
barDouble = sprintf('%s\n', repmat('=', 1, fieldWidth1 + fieldWidth2));
barSingle = sprintf('%s\n', repmat('-', 1, fieldWidth1 + fieldWidth2));
% re-organize the data 
c = c(2:end,[2 1 3 4 5]);
c = [num2cell(repmat(fieldWidth1, size(c,1), 1)), effect(2:end), c]';
% create the ANOVA table
s = sprintf('ANOVA table for %s\n', dvName);
s = [s barDouble];
s = [s sprintf('%-*s %4s %11s %14s %9s %9s\n', fieldWidth1, 'Effect', 'df', 'SS', 'MS', 'F', 'p')];
s = [s barSingle];
s = [s, sprintf('%-*s %4d %14.3f %14.3f %10.3f %10.4f\n', c{:})];
s = [s, barDouble];
end
