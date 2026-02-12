function [fn, failureFlag] = locateVertices(subj)

failureFlag = true; fn = []; 

load('/Users/diamondjm/Documents/Josh_work/Extraction/outcomesTable.mat','outcomesTable');

tableInd = strcmp(outcomesTable.ZaghloulNumber,subj);
if ~any(tableInd); return; end 


tableInd = outcomesTable.pNumber{tableInd};
if isempty(tableInd); return; end 

fn = '/Users/diamondjm/Documents/Josh_work/Masks/Price/';
fn = fullfile(fn,tableInd);

if exist(fn,'dir')==7; failureFlag = false; end
    
end
