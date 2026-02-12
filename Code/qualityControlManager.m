function [myStruct, wasModified] = qualityControlManager(myStruct,qualityControlThresh,varargin)

% Enforce quality control and presence of ROIs.


p = inputParser;
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
useAlternative = p.Results.useAlternative;

wasModified = false;

for kk = 1:length(myStruct)
    
    if isempty(myStruct(kk).subjID); continue; end
    % Really only necessary under the old paradigm. Recall that for NIH089
    % we started processing empty structs 
    
    [runQc,runRoi] = qualityControlChecker(myStruct(kk),qualityControlThresh);
    
    %%
    
    if runQc
        myStruct(kk) = qualityControlCycle(myStruct(kk),qualityControlThresh);
        wasModified = true; 
    end
    if runRoi
        % if isempty(myStruct(kk).sourceLoc.localizationResults); continue; end
        % We're now supporting localizing ROIs in empty structs. 
        myStruct(kk) = rectifyROI(myStruct(kk),'forceNewROI',true,'useAlternative',useAlternative);
        wasModified = true; 
    end
    
end