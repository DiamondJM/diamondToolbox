function plotDimRedVert2(mySingleton,varargin)



%% Preamble

assert(isscalar(mySingleton)); 
myBp = mySingleton.sourceLoc.braindata.myBp; 

currentVertices = mySingleton.sourceLoc.localizationResults(1,:);
locationsMaster = convertVerticesToLocations(mySingleton);

assert(~isempty(locationsMaster),'No localization to plot.'); 

p = inputParser;
addParameter(p,'colorVals',nan);
parse(p,varargin{:})
colorVals = p.Results.colorVals;

if isfield(mySingleton,'SeizureType'); eventType = 'seizure'; 
elseif isfield(mySingleton,'baselineNum'); eventType = 'baseline';
end


%% Process dimension reduction.


score = locationsMaster;
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

%% Plot 
% cla
hold on

% [resecPointLocs, failureFlag] = retrieveMaskValues(mySingleton.subjID);
[lhLoc,rhLoc] = resectedVerticesLocations(mySingleton);
resecPointLocs = [lhLoc; rhLoc];

resecPointLocs = resecPointLocs * coeff;
% resecPointLocs = (resecPointLocs - mu) * coeff;
resecPointLocs = double(resecPointLocs(:,[1 2]));
boundInds = boundary(resecPointLocs(:,1),resecPointLocs(:,2));
plot(resecPointLocs(boundInds,1),resecPointLocs(boundInds,2),'k-','linewidth',2);

resecPointLocs = [myBp.surfaces.pial_lh.vertices; myBp.surfaces.pial_rh.vertices];
resecPointLocs = resecPointLocs * coeff;
resecPointLocs = resecPointLocs(1:10:end,:);
resecPointLocs = double(resecPointLocs(:,[1 2]));
boundInds = boundary(resecPointLocs(:,1),resecPointLocs(:,2));
plot(resecPointLocs(boundInds,1),resecPointLocs(boundInds,2),'k-','linewidth',2);

for ii = 1:size(score,1)
    plot(score(ii,1),score(ii,2),'.', ...
        'color',colorVals(ii,:),...
        'MarkerSize',countVals(ii))
end

%% Formatting

[d1, d0] = coeffToOrientation(coeff);

axis equal
set(gca,'DataAspectRatio',[1 1 1])
box on

xticks(xlim); yticks(ylim);
xticklabels({sprintf('%s',d0{1,1}), sprintf('%s',d1{1,1})});
yticklabels({sprintf('%s',d0{2,1}), sprintf('%s',d1{2,1})});

