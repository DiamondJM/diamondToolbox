function sourceLoc = getSensors(sourceLoc,varargin)

p = inputParser;
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
useAlternative = p.Results.useAlternative; 


sensorDistance = sourceLoc.paramStruct.sensorDistance; 

geodesic = loadGeodesic(sourceLoc,'useAlternative',useAlternative); 
geodesicMaster = geodesic.geodesicDistances;
vertexNums = geodesic.vertexNums;

isLeft = geodesic.isLeftInds;

% distancesMaster = geodesicMaster(:,vertexNums);
% Unfortunately the above doesn't work anymore, now with Depths.
% Some elements of vertexNums are NaNs, corresponding to depth electrodes I
% couldn't localize. 
distancesMaster = nan(length(vertexNums)); 
isGood = ~isnan(vertexNums); 
distancesMaster(:,isGood) = geodesicMaster(:,vertexNums(isGood));

oppositeSideMask = xor(isLeft,isLeft');
distancesMaster(oppositeSideMask) = nan;
% distancesMaster(oppositeSideMask) = Inf; 
% warning('Contralateral sensors Inf.'); 

[ii,jj] = find(distancesMaster <= sensorDistance);
sensorsMaster = [ii jj];

sensorsMaster = sensorsMaster(ii < jj,:);

sourceLoc.sensors.sensorInds = sensorsMaster;
sourceLoc.sensors.distancesMaster = distancesMaster;
sourceLoc.sensors.sensorDistance = sensorDistance;

end