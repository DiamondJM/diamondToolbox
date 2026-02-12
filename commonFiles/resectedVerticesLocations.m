function [lhLoc,rhLoc,failureFlag] = resectedVerticesLocations(subj,varargin)

p = inputParser;
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
useAlternative = p.Results.useAlternative;

[~,subj] = rectifySubj(subj); 

[~, myBp] = retrieveBraindata(subj,'useAlternative',useAlternative);
% load('/Users/diamondjm/Documents/Josh_work/commonFiles/Legacy/Temp_BD/NIH023/bdBpFull.mat')

[rh,lh,failureFlag] = retrieveResectedVertex(subj);

lhVert = myBp.surfaces.pial_lh.vertices;
rhVert = myBp.surfaces.pial_rh.vertices;

lhLoc = lhVert(lh,:); 
rhLoc = rhVert(rh,:); 

end