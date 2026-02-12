function [goodOutcome, failureFlag,fullOutcome] = retrieveOutcome(subj)


[~,subj] = rectifySubj(subj); 
% failureFlag = true;

% We'll let NaN mean as yet unknown outcome. 
% 0, on the other hand, means outcome irretrievable, in other words, lost
% to follow-up or did not have resection. 

% Also, the function will work as follows. 
% outcome will just be a bool. 
% 1 means good outcome. 
% 0 means bad or no outcome. 
% And then failureFlag will be used to designate no outcome. 

% Temp 
load('/Users/diamondjm/Documents/Josh_work/Extraction/outcomesTable.mat','outcomesTable')


tableInd = strcmpi(outcomesTable.ZaghloulNumber,subj);
outcomeExists = outcomesTable.outcome; 

if any(tableInd)

    if isnan(outcomeExists(tableInd)) || outcomeExists(tableInd) == 0
        failureFlag = true; 
        goodOutcome = false; 
    else
        fullOutcome = outcomesTable.engelString(tableInd);
        % goodOutcome = strcmp(outcome,'1a');
        goodOutcome = outcomeExists(tableInd) == 1; % For 1a or 1b 
        
        % goodOutcome = outcomeExists(tableInd) <= 3; warning('Using Engel 3 as good.'); 
        
        % goodOutcome = outcomeExists(tableInd); 
        failureFlag = false; 
    end
else
    failureFlag = true; goodOutcome = 0; 
end

if failureFlag
    fprintf('Failed to find outcome for %s. \n',subj);
end

end
