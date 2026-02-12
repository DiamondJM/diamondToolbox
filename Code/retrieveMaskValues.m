function [maskValues, failureFlag] = retrieveMaskValues(ptString)

failureFlag = false;
failedKate = false; 


%%% First let's try using Kate's resection volumes.
load('/Users/diamondjm/Documents/Josh_work/Extraction/outcomesTable','outcomesTable');
tableInd = strcmpi(outcomesTable.ZaghloulNumber,ptString);
if ~any(tableInd); failedKate = true;
else
    cesString = outcomesTable.pNumber{tableInd};
    if ~isempty(cesString)
        cesDir = '/Users/diamondjm/Documents/Josh_work/Masks/Katie';
        % cesDir = '/Volumes/Shares/EEG/Users/IRTA/josh/rsxn_masks';
        cesContents = lsCell(cesDir);
        
        if ismember(cesString,cesContents); fileName = 'mvXfm';
        else; failedKate = true;
        end
    else
        failedKate = true;
    end
    
end

if failedKate
    
    load('/Users/diamondjm/Documents/Josh_work/Masks/Alison/conversionMaster','conversionMaster')
    tableInd = strcmpi(conversionMaster(:,2),ptString);
    
    if any(tableInd) && ~isequal(ptString,'NIH035') % Problematic patient 
        cesString = conversionMaster{tableInd,1};
        cesDir = '/Users/diamondjm/Documents/Josh_work/Masks/Alison/4_resec_mask/';        
        warning('Failed to find Kate''s resection mask for %s. We found Alison''s, which we''re using.',ptString)
        fileName = 'MaskValues';
        
    else
        fprintf('Failed to find mask for %s. Returning. \n',ptString)
        maskValues = [];
        failureFlag = true;
        return
    end
end

resecMask = fullfile(cesDir,cesString,fileName);
maskValues = dlmread(resecMask);

maskValues = maskValues(:,1:3);

if ~failureFlag
    warning('Mask found for %s, but from deprecated method.',ptString)
    % Function deprecated as of Price's new masks.
    % We should only get here if we're looking for resection info for a patient
    % Price couldn't process.

end

end