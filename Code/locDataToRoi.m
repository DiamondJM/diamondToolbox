function mySingleton = locDataToRoi(mySingleton,varargin)

%% Preamble
tic
assert(length(mySingleton) == 1,'Takes a single element of myStruct.')

p = inputParser;
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
useAlternative = p.Results.useAlternative; 


[myBd, ~] = retrieveBraindata(mySingleton.subjID,'useAlternative',useAlternative);

roiRadius = Inf;
if ~isinf(roiRadius); warning('ROI radius of %d',roiRadius); end

sumEmpty = 0;

%%

currentLocs = mySingleton.sourceLoc.localizationResults;
currentLocs(:,isnan(currentLocs(1,:))) = [];

if isempty(currentLocs)
    vertexMap = []; % Empty 
else
    
    currentVertices = currentLocs(1,:);
    
    uniqueVertex = unique(currentVertices);
    vertexMap = containers.Map('keyType','double','valueType','any');
    
    %%
    
    for vertIndex = uniqueVertex
        currentInds = find(currentVertices==vertIndex);
        numVertices = length(currentInds);
        
        isLeft = sign(vertIndex) == -1;
        
        if isLeft; surfString = 'lh'; currentSign = -1;
        else; surfString = 'rh'; currentSign = 1;
        end
        % Let's give left-sided vertices a negative sign.
        
        [roicUnique, ~, roiOutput] = myBd.vertex2ROI(abs(vertIndex),surfString,roiRadius);
        if isempty(roicUnique); sumEmpty = sumEmpty + numVertices; continue; end
        
        roicUnique = revertToVector(roicUnique);
        roicAll = roiOutput.ROIC_mesh_ndx;
        
        for roicIndex = roicUnique
            
            %% Upload, second-wise
            
            if isKey(vertexMap,roicIndex * currentSign)
                roiStruct = vertexMap(roicIndex * currentSign);
            else
                roiStruct = struct;
                roiStruct.count = 0;
                roiStruct.roi = roiOutput.vertex(roicAll == roicIndex);
                roiStruct.times = zeros(1,0);
            end
            
            roiStruct.count = roiStruct.count + numVertices;
            roiStruct.times = [roiStruct.times currentLocs(3,currentInds)];
            
            vertexMap(roicIndex * currentSign) = roiStruct;
            
        end
    end
    
    if sumEmpty; fprintf('%d unrecognized vertices found for %s.\n',sumEmpty,mySingleton.subjID); end
    % fprintf('ROI data computed in %.1f seconds. \n',toc)

end

roiResults = struct;
roiResults.vertexMap = vertexMap;
roiResults.timeWindow = 0;
roiResults.roiRadius = roiRadius;

if isfield(mySingleton.sourceLoc,'paramStruct') && isfield(mySingleton.sourceLoc.paramStruct,'qualityControlThresh')
    roiResults.qualityControlThresh = mySingleton.sourceLoc.paramStruct.qualityControlThresh;
end

mySingleton.sourceLoc.roiResults = roiResults;


