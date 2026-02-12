function vert = convertLocationsToVertices(myBp,locations)

vertL = myBp.surfaces.pial_lh.vertices; 
vertR = myBp.surfaces.pial_rh.vertices; 

stdVert = myBp.stdNumVert; 

vertAll = [vertL;vertR]; 

[~,vert] = min(pdist2(vertAll,locations)); 

isLeft = vert <= stdVert;
vert(isLeft) = vert(isLeft) * -1; 
vert(~isLeft) = vert(~isLeft) - stdVert; 

end