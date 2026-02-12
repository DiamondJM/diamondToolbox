function sourceLoc = getSensors(sourceLoc)

sensorDistance = sourceLoc.paramStruct.sensorDistance; 

geodesic = sourceLoc.geodesic;
vertexNums = geodesic.vertexNums;

isLeft = geodesic.isLeftInds;

distancesMaster = geodesic.geodesicDistances(:,vertexNums);
oppositeSideMask = xor(isLeft,isLeft');
distancesMaster(oppositeSideMask) = nan;
% Prevent getting sensor pairs on opposite sides of the brain. There's no
% way for us to define the geodesic distance between such pairs. 

[ii,jj] = find(distancesMaster < sensorDistance);
sensorsMaster = [ii jj];

sensorsMaster = sensorsMaster(ii < jj,:);

sourceLoc.sensors.sensorInds = sensorsMaster;
sourceLoc.sensors.distancesMaster = distancesMaster;
sourceLoc.sensors.sensorDistance = sensorDistance;

end