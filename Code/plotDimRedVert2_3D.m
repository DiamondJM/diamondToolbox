function plotDimRedVert2_3D(mySingleton)



%% Preamble

currentVertices = mySingleton.sourceLoc.localizationResults(1,:);
locationsMaster = convertVerticesToLocations(mySingleton);
% Need to make sure the above plays nicely with everything else. 


%% Process dimension reduction.


% coeff = eye(3); 
% mu = pcaInfo.mu;

score = locationsMaster; 
% score = (locationsMaster - mu) * coeff;

inRangeX = arrayfun(@(x) sum(abs(score(:,1) - score(x,1)) <= 3),1:size(score,1));
inRangeY = arrayfun(@(x) sum(abs(score(:,2) - score(x,2)) <= 3),1:size(score,1));
inRange = inRangeY .* inRangeX;


[a,b,c] = unique(currentVertices); 
score = score(b,:); 
locsCounts = histcounts(c,'binEdges',1:1:length(b) + 1);    
    
% The above should index into a, where each element is the number of
% appearances of a. 

% We can sort by size so we plot largest first 
[~,sortCount] = sort(locsCounts,'descend'); 
score = score(sortCount,:); 
locsCounts = locsCounts(sortCount); 
countVals = normalizeToBounds(locsCounts,[60 120]);  % Scaling

inRange = inRange(b);
inRange = inRange(sortCount);
cmap = jet;

colorVals = round(normalizeToBounds(inRange,[1 size(cmap,1)]));

% ColorVals may be passed in by the function
% plotEarlySourceDimRedWrapper.

% colorVals = cmap(colorVals,:);

%% Plot 

sourceLoc = mySingleton.sourceLoc;

currentEl = -90;

%%%% 
% For seizure source info
%%%%

[myBd, myBp] = retrieveBraindata(mySingleton.subjID);


geodesic = loadGeodesic(sourceLoc); 
isLeftInds = geodesic.isLeftInds;
useLeft = logical(round(sum(isLeftInds) / length(isLeftInds)));


myBd.ezplot(myBp,gca);
% plotResectionSurf(subj)
if useLeft; view(-90,currentEl);
else; view(90,currentEl);
end
myBp.camlights(5);

hold on; 
for ii = 1:size(score,1)
    plot3(score(ii,1),score(ii,2),score(ii,3),'.', ...
        'color',cmap(colorVals(ii),:), ...
        'MarkerSize',countVals(ii))
end

end

