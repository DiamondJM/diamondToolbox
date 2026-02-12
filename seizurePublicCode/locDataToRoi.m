function mySingleton = locDataToRoi(mySingleton,timeWindow)

%% Preamble 
tic

sourceLoc = mySingleton.sourceLoc;
lengthTs = sourceLoc.lengthTs; 

myBd = sourceLoc.braindata.myBd; 

localizationResults = sourceLoc.localizationResults;

if nargin == 1; timeWindow = 1000; 
else; if ~timeWindow; timeWindow = lengthTs; end 
end

stepsBuffer = buffer(1:lengthTs,timeWindow);
vMapAll = cell(size(stepsBuffer,2),1);

maxValAll = zeros(1,size(stepsBuffer,2)); 

roiRadius = Inf;

sumEmpty = 0;
parfor jj = 1:size(stepsBuffer,2)
    
    maxVal = 0; 

    currentTime = stepsBuffer(:,jj);
    currentTime = currentTime(logical(currentTime)); 
    
    currentLocs = localizationResults(:,ismember(localizationResults(3,:),currentTime));
    
    currentVertices = currentLocs(1,:);
    
    if isempty(currentVertices); continue; end    
    uniqueVertex = unique(currentVertices);
    
    vertexMap = containers.Map('keyType','double','valueType','any');
    
    for vertIndex = uniqueVertex
        currentInds = find(ismember(currentVertices,vertIndex));
        numVertices = length(currentInds);
        
        isLeft = sign(vertIndex) == -1; 
        
        if isLeft; surfString = 'lh'; currentSign = -1;
        else; surfString = 'rh'; currentSign = 1;
        end
        % Let's give left-sided vertices a negative sign. 
        
        [roicUnique, ~, roiOutput] = myBd.vertex2ROI(abs(vertIndex),surfString,roiRadius);
        
        if isempty(roicUnique); sumEmpty = sumEmpty + numVertices; continue; end
        roicUnique = roicUnique';
        roiAll = roiOutput.ROIC_mesh_ndx;
        
        for roicIndex = roicUnique
            
            if isKey(vertexMap,roicIndex * currentSign)
                roiStruct = vertexMap(roicIndex * currentSign);
                roiStruct.count = roiStruct.count + numVertices;
            else
                roiStruct = struct;
                roiStruct.count = numVertices;
                roiStruct.roi = roiOutput.vertex(roiAll == roicIndex);
            end
                 
            maxVal = max(maxVal,roiStruct.count);
            vertexMap(roicIndex * currentSign) = roiStruct;
        end        
        
    end
        
    vMapAll{jj} = vertexMap;
    
    maxValAll(jj) = maxVal; 
    
end

fprintf('%d unrecognized vertices found. \n',sumEmpty); 
fprintf('ROI data computed in %.1f seconds. \n',toc)

roiResults = struct;
roiResults.vMapAll = vMapAll; 
roiResults.timeWindow = timeWindow; 
roiResults.maxValAll = maxValAll; 

mySingleton.sourceLoc.roiResults = roiResults; 

