function [lhLoc,rhLoc] = resectedVerticesLocations(mySingleton)


% myBd = mySingleton.sourceLoc.braindata.myBd; 
myBp = mySingleton.sourceLoc.braindata.myBp; 

lh = mySingleton.sourceLoc.resectedVertices.lh; 
rh = mySingleton.sourceLoc.resectedVertices.rh; 

lhVert = myBp.surfaces.pial_lh.vertices;
rhVert = myBp.surfaces.pial_rh.vertices;

lhLoc = lhVert(lh,:); 
rhLoc = rhVert(rh,:); 

end