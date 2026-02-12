function mySingleton = rectifyROI(mySingleton,timeWindow)

%% Preamble

forceNewROI = false; 
% warning('Forcing new ROI.'); 

lengthTs = mySingleton.sourceLoc.lengthTs; 

if nargin == 1
    timeWindow = 1000; 
    % Default timeWindow
else
    if isfield(mySingleton.sourceLoc,'roiResults')
        if ~forceNewROI
            roiResults = mySingleton.sourceLoc.roiResults;
            if roiResults.timeWindow == lengthTs; return; end
            if roiResults.timeWindow == timeWindow; return; end
        end
    end
end


mySingleton = locDataToRoi(mySingleton,timeWindow);
% mySingleton = locDataToRoi_MaxCum(mySingleton,timeWindow);
