function plotResectedVertices(subj,varargin)
%% Preamble 

p = inputParser;
addParameter(p,'newPlot',true);
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
newPlot = p.Results.newPlot;
useAlternative = p.Results.useAlternative;

[~,subj] = rectifySubj(subj); 

%% Plotting

[myBd, myBp] = retrieveBraindata(subj,'useAlternative',useAlternative);
% load('/Users/diamondjm/Documents/Josh_work/commonFiles/Legacy/Temp_BD/NIH023/bdBpFull.mat')
if newPlot; myBd.ezplot(myBp,gca); end

[lhLoc,rhLoc] = resectedVerticesLocations(subj,'useAlternative',useAlternative);

plot3(rhLoc(:,1),rhLoc(:,2),rhLoc(:,3),'k.','markersize',20);
% plot3(lhLoc(:,1),lhLoc(:,2),lhLoc(:,3),'b.','markersize',10);
plot3(lhLoc(:,1),lhLoc(:,2),lhLoc(:,3),'k.','markersize',20);

% if useLeft; view(-90,currentEl);
% else; view(90,currentEl);
% end
myBp.camlights(5);

