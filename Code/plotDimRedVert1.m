function  pcaInfo = plotDimRedVert1(mySingleton,varargin)
%% Grab and process dimension reduction.

locationsMaster = convertVerticesToLocations(mySingleton);
if size(locationsMaster,1) == 0
    warning('No points to plot, patient %s, epoch %s.',mySingleton.subjID,mySingleton.eegfile(end-10:end)); 
    return; 
end

p = inputParser;
addParameter(p,'colorVals',repmat([0 0 0],size(mySingleton.sourceLoc.localizationResults,2),1));
addParameter(p,'pcaInfo',generatePcaInfo(locationsMaster));
addParameter(p,'satModifier',nan); 

parse(p,varargin{:})
colorVals = p.Results.colorVals;
pcaInfo = p.Results.pcaInfo; 
satModifier = p.Results.satModifier; 


coeff = pcaInfo.coeff; mu = pcaInfo.mu; 
% At some point I should standardize dimred1 and 2, so that both have the
% same handling of mu, coeff, pcaInfo, et cetera. 
% Why do we use mu here but in in dimension 2? 

score = (locationsMaster - mu) * coeff;
timeVals = mySingleton.sourceLoc.localizationResults(3,:) / 1000; 


if ~ismember('satModifier',p.UsingDefaults)
    
    [~,inds] = sort(satModifier,'ascend');
    colorVals = colorVals(inds,:);
    score = score(inds,:);
    timeVals = timeVals(inds);
end
%% Plot 

[lh,rh, failureFlag] = retrieveResectedVertex_orRoi(mySingleton.subjID);
resecPointLocs = [find(rh); -find(lh)];
resecPointLocs = convertVerticesToLocations(mySingleton,revertToVector(resecPointLocs)); 
if ~failureFlag
    resecPointLocs = (resecPointLocs - mu) * coeff;
    [a,b] = bounds(resecPointLocs(:,1));
else
    [a,b] = deal(0); 
end

% Create our figure  
% cla
hold on

% colorVals in this setting is called by evalSpkClusters.
for ii = 1:size(score,1)
    plot(timeVals(ii),score(ii,1),'.','MarkerSize',20,'color',colorVals(ii,:));
end

%% Formatting

[d1, d0] = coeffToOrientation(coeff);

tempXlim = xlim;

plot([tempXlim(2) tempXlim(2)],[a b],'k','linewidth',3)

yticks(ylim);
yticklabels([d0(1,1) d1(1,1)]);
% tempYticks = yticklabels; 
% tempYticks(1) = d0(1,1); 
% tempYticks(end) = d1(1,1); 
% yticklabels(tempYticks); 

tempYlim = ylim;
if tempYlim(2) > 0 && tempYlim(1) < 0; plot([0 0],tempYlim,'k.'); end


%% Pack up 

pcaInfo.mu = mu;
pcaInfo.coeff = coeff; 

end

