function [tagNames,chansToTag,colorsAll,colorList] = groupByHardware(sourceLoc)

chansMaster = generateChansMasterList(sourceLoc.subj);

chansToTag = cellfun(@(x) isstrprop(x,'digit'),chansMaster,'UniformOutput',0);
chansToTag = cellfun(@(x,y) x(~y),chansMaster,chansToTag,'UniformOutput',0);


tagNames = unique(chansToTag); 

[~,chansToTag] = ismember(chansToTag,tagNames);
 

%% Colors 

colorList = parula(length(tagNames));
% colorsAll = colorList(chansToTag,:);
% Doesn't work if some channels are missing. 
colorsAll = nan(length(chansToTag),3); 

for ii = 1:length(tagNames)
    currentInds = chansToTag == ii; 
    % currentColors = colorsAll(currentInds,:); 
    currentColors = repmat(colorList(ii,:),sum(currentInds),1); 
    
    currentColors = rgb2hsv(currentColors); 
    
    satModifier = linspace(currentColors(1,2),.2,sum(currentInds))'; 
    brightModifier = linspace(currentColors(1,3),.8, sum(currentInds))'; 
    
    currentColors(:,2) = satModifier; 
    currentColors(:,3) = brightModifier; 
    
    currentColors = hsv2rgb(currentColors); 
    
    currentText = cellfun(@(x) str2double(x),erase(chansMaster(currentInds),tagNames(ii)));
    [~,ind] = sort(currentText,'ascend'); 
    currentColors = currentColors(ind,:); 
    
    colorsAll(currentInds,:) = currentColors; 
end

%% Rectify to chanNames 

chanNames = sourceLoc.chanNames; 

[~,inds] = ismember(chanNames,chansMaster);
assert(all(logical(inds)),'chanNames passed which we don''t have in our records.'); 

chansToTag = chansToTag(inds); colorsAll = colorsAll(inds,:); 


end