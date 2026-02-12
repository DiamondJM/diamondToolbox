function myStruct = plotSurfFun(myStruct,timeWindow)


%% Preamble 

sourceLoc = myStruct.sourceLoc;

currentEl = -90;

%%%% 
% For seizure source info
%%%%

myBd = sourceLoc.braindata.myBd; 
myBp = sourceLoc.braindata.myBp; 

geodesic = myStruct.sourceLoc.geodesic;
isLeftInds = geodesic.isLeftInds;
useLeft = logical(round(sum(isLeftInds) / length(isLeftInds)));

if nargin == 1; timeWindow = 1000; end
% If timeWindow = 0, we can plot localization results over the course of
% the entire seizure, on a single brain plot. 
% Otherwise, we can create a source localization video, where results are
% plotted over windows of length timeWindow (samples). 
% The default timeWindow is of length 1 second. 

%% Brain plot 

%%%%%
% Basic setup
%%%%%

colormap(jet);
% myBd.ezplot(myBp,gca); % If we would not like to include the resection territory 
plotResectionSurf(myStruct) % If we would like to include the resection territory 
if useLeft; view(-90,currentEl);
else; view(90,currentEl);
end
myBp.camlights(5);
ax = gca;

%%  Grab and set up ROI results

myStruct = rectifyROI(myStruct,timeWindow);
roiResults = myStruct.sourceLoc.roiResults; 
vMapAll = roiResults.vMapAll; 
% clim = roiResults.clim; 


%% Find max color 

% The coloring is subject to the amount of ROI hits and so is difficult to
% predict before mapping the data to ROIs. Therefore first we determine max
% color, and then we plot. 

% This section is only needed for timeWindow ~= 0. 

maxColor = 0; 
for jj = 1:length(vMapAll)
    if isempty(vMapAll{jj}); continue; end
        
    vertexMap = vMapAll{jj};
    
    mapKeys = cell2mat(vertexMap.keys);
    VPerROI = cell(size(mapKeys));
    valPerROI = zeros(size(mapKeys));
    
    for ii = 1:length(mapKeys)
        roiStruct = vertexMap(mapKeys(ii));
        VPerROI{ii} = roiStruct.roi;
        valPerROI(ii) = roiStruct.count;
    end
    isLeft = sign(mapKeys) == -1;
    
    if all(isLeft)
        currentMax = findMaxColor(myBp, VPerROI, valPerROI,'surf','lh');
    elseif all(~isLeft)
        currentMax = findMaxColor(myBp, VPerROI, valPerROI,'surf','rh');
    else
        [isLeft,sortInds] = sort(isLeft,'descend');
        VPerROI = VPerROI(sortInds);
        valPerROI = valPerROI(sortInds);
        rh_begin = find(~isLeft,1);
        
        currentMax = findMaxColor(myBp, VPerROI, valPerROI,'rh_begin',rh_begin);
        
    end
    maxColor = max(maxColor,currentMax);    
end
fprintf('Max color is %.f. \n',maxColor); 

clim = [0 maxColor];

%% Plot what we have 
    
lastPlot = false;

for jj = 1:length(vMapAll) 

    if ~isempty(vMapAll{jj})
        
        vertexMap = vMapAll{jj};
        
        mapKeys = cell2mat(vertexMap.keys);
        VPerROI = cell(size(mapKeys));
        valPerROI = zeros(size(mapKeys)); 
        
        for ii = 1:length(mapKeys)
            roiStruct = vertexMap(mapKeys(ii));
            VPerROI{ii} = roiStruct.roi; 
            valPerROI(ii) = roiStruct.count; 
        end
        
        
        isLeft = sign(mapKeys) == -1; 
        
        axes(ax);
        if lastPlot; myBp.clearRegions(); end
        
        if all(isLeft)
            myBp.plotRegionsData(VPerROI, valPerROI,'surf','lh','clim',clim);
        elseif all(~isLeft)
            myBp.plotRegionsData(VPerROI, valPerROI,'surf','rh','clim',clim);
        else 
            [isLeft,sortInds] = sort(isLeft,'descend'); 
            VPerROI = VPerROI(sortInds); 
            valPerROI = valPerROI(sortInds); 
            rh_begin = find(~isLeft,1);
            
            myBp.plotRegionsData(VPerROI, valPerROI,'rh_begin',rh_begin,'clim',clim);
        end
        fprintf('%d ROICs plotted; %d total ROI hits. \n',length(mapKeys),sum(valPerROI));   
        lastPlot = true; 
    elseif lastPlot 
        axes(ax);
        myBp.clearRegions();
        lastPlot = false; 
    end
    drawnow    
end


end