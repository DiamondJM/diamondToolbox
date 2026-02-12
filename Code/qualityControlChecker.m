function [runQc,runRoi] = qualityControlChecker(mySingleton,qualityControlThresh)

assert(length(mySingleton)==1); 

runQc = false;
runRoi = false;

%%

if isfield(mySingleton.sourceLoc.paramStruct,'qualityControlThresh') ...
        && mySingleton.sourceLoc.paramStruct.qualityControlThresh==qualityControlThresh
    % Do nothing for quality control
    
    % Check ROI results
    if isfield(mySingleton.sourceLoc,'roiResults') ...
            && isfield(mySingleton.sourceLoc.roiResults,'qualityControlThresh') ...
            && mySingleton.sourceLoc.roiResults.qualityControlThresh == qualityControlThresh
        % Do nothing, we're good
    else
        runRoi = true;
    end
else
    runQc = true;
    runRoi = true;
end

end