function myStruct = qualityControlCycle(myStruct,qualityControlThresh)

if nargin == 1; qualityControlThresh = Inf; end

useCycle = true; 

for kk = 1:length(myStruct)
    
    assert(isfield(myStruct(kk),'sourceLoc'));
    if ~isfield(myStruct(kk).sourceLoc.paramStruct,'locResultsOrig'); useCycle = false; end
    % If it is a field, useCycle remains true
    % So long as useCycle is true, it's impossible to overwrite the field 
    
    if useCycle
        % Cycle
        if myStruct(kk).sourceLoc.paramStruct.qualityControlThresh == qualityControlThresh; continue; end
        localizationResults = myStruct(kk).sourceLoc.paramStruct.locResultsOrig;
    else
        localizationResults = myStruct(kk).sourceLoc.localizationResults;
        myStruct(kk).sourceLoc.paramStruct.locResultsOrig = localizationResults;
        % This is the key step that overwrites the field, and can only be
        % accessed if useCycle is off. 
        % Use cycle is true by default. 
        
        % Furthermore, as of 10/9/22, we've started adding locResultsOrig
        % and qualityControlThresh in localizationFunction. 
        % So we should never get here, in the future. 
        warning('I got here.'); 
    end
    
    %% Process 
    
    if ~isinf(qualityControlThresh)
        badInds = localizationResults(2,:) > qualityControlThresh;
        fprintf('%.2f%% of points removed in quality control, for %d total. \n',sum(badInds) / length(badInds) * 100,sum(badInds));
        localizationResults(:,badInds) = [];
    end

    %% Pack 
    
    myStruct(kk).sourceLoc.localizationResults = localizationResults;
    myStruct(kk).sourceLoc.paramStruct.qualityControlThresh = qualityControlThresh;
end

end