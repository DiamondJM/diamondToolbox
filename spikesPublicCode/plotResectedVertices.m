function plotResectedVertices(mySingleton)


%% Plotting

myBd = mySingleton.sourceLoc.braindata.myBd; 
myBp = mySingleton.sourceLoc.braindata.myBp; 

myBd.ezplot(myBp,gca);

[lhLoc,rhLoc] = resectedVerticesLocations(mySingleton);

plot3(rhLoc(:,1),rhLoc(:,2),rhLoc(:,3),'k.','markersize',10);
plot3(lhLoc(:,1),lhLoc(:,2),lhLoc(:,3),'b.','markersize',10);

% if useLeft; view(-90,currentEl);
% else; view(90,currentEl);
% end
myBp.camlights(5);

