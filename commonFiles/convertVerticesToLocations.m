function locationsMaster = convertVerticesToLocations(mySingleton,localizationResults,myBp,varargin) 
%% Preamble 


p = inputParser;
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
useAlternative = p.Results.useAlternative;


% assert(length(mySingleton) == 1,'Takes a single element of myStruct, %s.',mySingleton(1).subjID)
assert(length(mySingleton) == 1);

if nargin <= 2 || isempty(myBp); [~, myBp] = retrieveBraindata(mySingleton.subjID,'useAlternative',useAlternative);end

lVertices = myBp.surfaces.pial_lh.vertices;
rVertices = myBp.surfaces.pial_rh.vertices;

%% Let's distribute our data into a master map. 

if nargin == 1 || isempty(localizationResults); localizationResults = mySingleton.sourceLoc.localizationResults; end 

locationsMaster = nan(size(localizationResults,2),3);

for ii = 1:size(localizationResults,2)
    isLeft = sign(localizationResults(1,ii)) == -1;
    if isLeft; currentVert = lVertices;
    else; currentVert = rVertices;
    end
    
    if isnan(localizationResults(1,ii)); continue; end
    
    locationsMaster(ii,:) = currentVert(abs(localizationResults(1,ii)),:);
end
        
end