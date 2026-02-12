function mySingleton = rectifyROI(mySingleton,varargin)

%% Preamble

p = inputParser;
% addParameter(p,'timeWindow',0);
addParameter(p,'forceNewROI',false);
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
% timeWindow = p.Results.timeWindow; 
forceNewROI = p.Results.forceNewROI;
useAlternative = p.Results.useAlternative; 

% lengthTs = mySingleton.sourceLoc.lengthTs; 
% lengthTs = max(mySingleton.sourceLoc.localizationResults(3,:)); 

%% Onto the function 

if ~isfield(mySingleton.sourceLoc,'roiResults') || forceNewROI
    mySingleton = locDataToRoi(mySingleton,'useAlternative',useAlternative); 
end

%% Old

% if isfield(mySingleton.sourceLoc,'roiResults') && ~forceNewROI
%     roiResults = mySingleton.sourceLoc.roiResults;
%     if roiResults.timeWindow == lengthTs
%         if timeWindow ~= 0 && timeWindow ~= lengthTs
%             forceNewROI = true;
%             % Prior time window was zero or lengthTs, resulting in lengthTs timeWindow.
%             % However, current timeWindow is neither.
%         end            
%     else
%         if roiResults.timeWindow ~= timeWindow; forceNewROI = true; end
%     end
%     
% else; forceNewROI = true; 
% end
% 
% if forceNewROI; mySingleton = locDataToRoi(mySingleton,'timeWindow',timeWindow); end

