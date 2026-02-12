function myBp = plotSensorsByHardware(subj,varargin)

%% Preamble 

sourceLoc.subj = subj; 

p = inputParser;
addParameter(p,'interestLeads',generateChansMasterList(subj));
addParameter(p,'newPlot',true); 
addParameter(p,'useAlternative',false); 
addParameter(p,'usePial',false);
parse(p,varargin{:})
interestLeads = p.Results.interestLeads; 
newPlot = p.Results.newPlot; 
useAlternative = p.Results.useAlternative;
usePial = p.Results.usePial; 

[myBd, myBp] = retrieveBraindata(subj,'useAlternative',useAlternative);


if newPlot; myBd.ezplot(myBp,gca); end

sourceLoc.chanNames = interestLeads; 
geodesic = loadGeodesic(sourceLoc,'useAlternative',useAlternative);
if useAlternative || usePial; leadLocs = geodesic.leadLocationsPial; % Only this is available for alternative surfaces
else; leadLocs = geodesic.leadLocations;
end

%% Determine channels list 

[tagNames,chansToTag,colorsAll,colorList] = groupByHardware(sourceLoc);


%% Take care of plotting 

if round(nnz(geodesic.isLeftInds)); xVal = -20; % Left-sided
else; xVal = 20; 
end

for ii = revertToVector(unique(chansToTag))

    modVal = [xVal,0,0];
    currentInds = find(chansToTag == ii); 
    colorsLocal = colorsAll(currentInds,:); 
    
    for jj = 1:length(currentInds)
        plot3(leadLocs(currentInds(jj),1),leadLocs(currentInds(jj),2),leadLocs(currentInds(jj),3),'.','MarkerSize',80,'color',colorsLocal(jj,:));
    end
    
    currentText = erase(interestLeads(currentInds),tagNames(ii));
    text(leadLocs(currentInds,1),leadLocs(currentInds,2),leadLocs(currentInds,3),currentText,'FontSize',18);
    
    [~,antInd] = max(leadLocs(currentInds,2));
    antInd = currentInds(antInd); 
    
    text(leadLocs(antInd,1) + modVal(1),leadLocs(antInd,2) + modVal(2),leadLocs(antInd,3) + modVal(3),tagNames(ii),'FontSize',42,'color',colorList(ii,:));
end

% for ii = 1:size(sensors,1)
%     sensorLine = leadLocs(sensors(ii,:),:);
%     
%     plot3(sensorLine(:,1),sensorLine(:,2),sensorLine(:,3),'-','color',[1 1 1],'LineWidth',3);
% end
if newPlot; myBp.setOpacity(0.75); end


end