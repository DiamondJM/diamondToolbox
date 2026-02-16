function pcaInfo = plotDimRedVert2(mySingleton,varargin)

%% Preamble

assert(isscalar(mySingleton)); 
% Need to make sure the above plays nicely with everything else. 

p = inputParser;
addParameter(p,'colorVals',nan);
addParameter(p,'pcaInfo',nan);
addParameter(p,'saving',false);
addParameter(p,'useAlternative',false);

% addParameter(p,'eventType','baseline')
% chanNames should now be passed as pair, value. 
addParameter(p,'useEye',true);
addParameter(p,'plotBoundaries',true); 

parse(p,varargin{:})
colorVals = p.Results.colorVals;
pcaInfo = p.Results.pcaInfo; 
useEye = p.Results.useEye;
useAlternative = p.Results.useAlternative;

[~, myBp] = retrieveBraindata(mySingleton.subjID,'useAlternative',useAlternative);

currentVertices = mySingleton.sourceLoc.localizationResults(1,:);
locationsMaster = convertVerticesToLocations(mySingleton,[],[],'useAlternative',useAlternative);

if isempty(locationsMaster)
    pcaInfo = struct; 
    warning('No localization results for %s, %s.',mySingleton.subjID,mySingleton.eegfile(end-10:end))
    return
end 


if ismember('pcaInfo',p.UsingDefaults); pcaInfo = generatePcaInfo(locationsMaster,'useEye',useEye); end 
% Sentinel value. Couldn't really find a better way to do this
% The problem with just putting the function call in the default value is
% that I don't really want to call the function every time I'm just
% passing something in anyway. 
saving = p.Results.saving; 
% eventType = p.Results.eventType; 
plotBoundaries = p.Results.plotBoundaries; 

eventType = mySingleton.task; 


%% Process dimension reduction.


coeff = pcaInfo.coeff;
% coeff = eye(3); 
% mu = pcaInfo.mu;

score = locationsMaster * coeff;
% score = (locationsMaster - mu) * coeff;
[a,b,c] = unique(currentVertices,'stable'); 
score = score(b,:); 
    
% The above should index into a, where each element is the number of
% appearances of a. 

% We can sort by size so we plot largest first 

% locsCounts = histcounts(c,'binEdges',1:1:length(b) + 1);    
% [~,sortCount] = sort(locsCounts,'ascend'); 
sortCount = randperm(size(score,1)); 
% sortCount = 1:size(score,1); 
% sortCount = size(score,1):-1:1; 
score = score(sortCount,:); 
% locsCounts = locsCounts(sortCount); 
% countVals = normalizeToBounds(locsCounts,[8 100]);  % Scaling

inRange = squareform(pdist(locationsMaster));
currentRad = 2; 
inRange = sum(inRange <= currentRad);

% maxDensity = max(inRange) / (4/3 * pi * currentRad ^ 3);
maxDensity = max(inRange);
fprintf('Max density is %.2f.\n',maxDensity); 

inRange = inRange(b);
inRange = inRange(sortCount);

% countVals = normalizeToBounds(inRange,[15 50],[1 108]); warning('Scaling.'); 
countVals = normalizeToBounds(inRange,[15 50]); 

if ismember('colorVals',p.UsingDefaults)
    
    if isequal(eventType,'baseline')
        cmap = turbo;
        % inRange = squareform(pdist(locationsMaster));
        % inRange = sum(inRange <= 2);
        
        % maxDensity = max(inRange) / (4/3 * pi * (2/2)^3);
        
        
        % inRange = inRange(b);
        % inRange = inRange(sortCount);
        
        colorVals = round(normalizeToBounds(inRange,[1 size(cmap,1)]));
        
        % ColorVals may be passed in by the function
        % plotEarlySourceDimRedWrapper.
        
        colorVals = cmap(colorVals,:);
    elseif isequal(eventType,'seizure')
        cmap = parula; 
        
        colorVals = zeros(size(a)); 
        
        for jj = 1:length(a)
            colorVals(jj) = mean(mySingleton.sourceLoc.localizationResults(3,mySingleton.sourceLoc.localizationResults(1,:)==a(jj))); 
        end
        colorVals = round(normalizeToBounds(colorVals,[1 size(cmap,1)]));
        colorVals = colorVals(sortCount); 
        
        
        % Below is optional I believe, the idea was to put the early points
        % on top 
        % [~,inds] = sort(colorVals,'descend'); 
        [~,inds] = sort(colorVals,'ascend'); 
        colorVals = colorVals(inds); 
        countVals = countVals(inds); 
        score = score(inds,:); 
        
        colorVals = cmap(colorVals,:);

    end
        
end



% In the spikes directory, there is no dimRedVert2. 
% That seems to be deprecated, as we're doing the size handling here. 
% Not necessarily bad, in fact probably better.
% But I don't explicitly recall why we dumped dimRedVert2. 
% But do note that this doesn't have the similar requirement to that in the
% seizure directory, which is that vertices with multiple localization
% events have to be colored according to the average time of localization. 

%% Plot 
% cla
hold on

if plotBoundaries
    % [resecPointLocs, failureFlag] = retrieveMaskValues(mySingleton.subjID);
    [lhLoc,rhLoc,failureFlag] = resectedVerticesLocations(mySingleton.subjID,'useAlternative',useAlternative);
    resecPointLocs = [lhLoc; rhLoc];
    
    if ~failureFlag
        resecPointLocs = resecPointLocs * coeff;
        % resecPointLocs = (resecPointLocs - mu) * coeff;
        resecPointLocs = double(resecPointLocs(:,[1 2]));
        boundInds = boundary(resecPointLocs(:,1),resecPointLocs(:,2));
        plot(resecPointLocs(boundInds,1),resecPointLocs(boundInds,2),'k-','linewidth',2);
    end
    
    resecPointLocs = [myBp.surfaces.pial_lh.vertices; myBp.surfaces.pial_rh.vertices];
    resecPointLocs = resecPointLocs * coeff;
    resecPointLocs = resecPointLocs(1:10:end,:);
    resecPointLocs = double(resecPointLocs(:,[1 2]));
    boundInds = boundary(resecPointLocs(:,1),resecPointLocs(:,2));
    plot(resecPointLocs(boundInds,1),resecPointLocs(boundInds,2),'k-','linewidth',2);
end

for ii = 1:size(score,1)
    plot(score(ii,1),score(ii,2),'.', ...
        'color',colorVals(ii,:),...
        'MarkerSize',countVals(ii))
end

%% Formatting

[d1, d0] = coeffToOrientation(coeff);

axis equal
set(gca,'DataAspectRatio',[1 1 1])
% set(gca,'PlotBoxAspectRatio',[3 4 4])
box on


xticks(xlim); yticks(ylim);
% xticklabels([d0(1,1) d1(1,1)]); 
% yticklabels([d0(2,1) d1(2,1)]);
xticklabels({sprintf('%s',d0{1,1}), sprintf('%s',d1{1,1})});
yticklabels({sprintf('%s',d0{2,1}), sprintf('%s',d1{2,1})});


% tempYlim = ylim;
% if tempYlim(2) > 0 && tempYlim(1) < 0; plot([0 0],tempYlim,'k.'); end

%% Wrap up 

pcaInfo.score = score; 
if saving
    
    fn = fullfile('Figures','dimRed'); 
    if isequal(eventType,'seizure'); fn = fullfile(fn,mySingleton.subjID); end
    
    if ~isfolder(fn); mkdir(fn); end 
    
    if isequal(eventType,'baseline'); saveas(gcf,fullfile(fn,sprintf('%s.svg',mySingleton.subjID)));  
    else; saveas(gcf,fullfile(fn,sprintf('%s.svg',sessionToIdentifier(mySingleton)))); 
    end

end

