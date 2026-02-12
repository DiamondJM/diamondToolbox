function locationsMaster = convertVerticesToLocations(mySingleton) 
%% Preamble 
 
assert(length(mySingleton) == 1,'Takes a single element of myStruct.')

sourceLoc = mySingleton.sourceLoc;
myBp = sourceLoc.braindata.myBp;

lVertices = myBp.surfaces.pial_lh.vertices;
rVertices = myBp.surfaces.pial_rh.vertices;

%% Let's distribute our data into a master map. 


localizationResults = mySingleton.sourceLoc.localizationResults;

locationsMaster = zeros(size(localizationResults,2),3);

for ii = 1:size(localizationResults,2)
    isLeft = sign(localizationResults(1,ii)) == -1;
    if isLeft; currentVert = lVertices;
    else; currentVert = rVertices;
    end
    
    locationsMaster(ii,:) = currentVert(abs(localizationResults(1,ii)),:);
end
        
end