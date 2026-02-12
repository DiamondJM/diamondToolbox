function geodesic = buildGeodesicDistances_master(subj,varargin)

%% Preamble 

p = inputParser;
addParameter(p,'useUtah',false);
addParameter(p,'useAlternative',false);
addParameter(p,'saving',true);

parse(p,varargin{:})
useUtah = p.Results.useUtah;
useAlternative = p.Results.useAlternative; 
saving = p.Results.saving; 

tic; 
[~,subj] = rectifySubj(subj);

[myBd, myBp, useAlternative] = retrieveBraindata(subj,'useAlternative',useAlternative);
% Pull back useAlternative; inside out methodology 

jacksheet = getJackTable(subj,pullServerDirectory);

fprintf('Building geodesic distances for %s. \n',subj); 
    
%% Go on to build the geodesic distances 

surfL = myBp.surfaces.pial_lh;
surfR = myBp.surfaces.pial_rh;
numVert = myBd.stdNumVertices;

if useAlternative
    
    % Start by loading prior; this must exist 
    sourceLoc.subj = subj;
    geodesic = loadGeodesic(sourceLoc,'unModified',true,'useAlternative',false,'useUtah',useUtah);
    
    isLeftInds = geodesic.isLeftInds;
    
    leadNames = geodesic.leadNames;
    numLeads = length(leadNames);
    distanceThresh = nan;
    distanceFull = nan(size(geodesic.distanceFull));
    
    
    if useUtah
        leadLocations = getUtahLocation(subj,'useAlternative',true);
        leadToVertex = convertLocationsToVertices(myBp,leadLocations); 
        
    else
        % The below all actually works JUST FINE with utah, too. But the
        % above, while less programatically elegant, is slightly more
        % precise. Because we have the exact utah location under the
        % alternative situation (because we picked it out).
        leadToVertex = geodesic.vertexNums;
        
        leadLocations = nan(size(geodesic.leadLocationsPial));
        % This is the slot for dural, but we don't have access to that information because it comes from tal.
        
    end

    
else
    
    if useUtah
        
        [~,t] = getUtahLocation(subj); 
    else
        
        leadsDir = fullfile(pullServerDirectory,subj,'/tal/leads.csv');
        assert(exist(leadsDir,'file'));
        t = readtable(leadsDir);
        
    end
    
    leadLocations = table2array(t(:,{'x','y','z'}));
    leadNames = t.chanName;
    numLeads = length(leadNames);
    
    %%%%%%%%%%%%%%
    % Which hemi?
    %%%%%%%%%%%%%%
    
    if useUtah
        
        vert = convertLocationsToVertices(myBp,leadLocations);
        isLeftInds = sign(vert)==-1;
        
    else
        [lLeads, rLeads] = myBd.chans2hem(leadNames, jacksheet);
        isLeftInds = ismember(leadNames,lLeads);
    end
    
    distanceThresh = 5; % mm
    
    leadToVertex = nan(1,numLeads);
    distanceFull = nan(1,numLeads);
    
    % Left
    leadDistance = pdist2(surfL.vertices,leadLocations(isLeftInds,:));
    [valsLeft,closestVert] = min(leadDistance);
    closestVert(valsLeft >= distanceThresh) = nan;
    leadToVertex(isLeftInds) = closestVert;
    distanceFull(isLeftInds) = valsLeft;
    
    % Right
    leadDistance = pdist2(surfR.vertices,leadLocations(~isLeftInds,:));
    [valsRight,closestVert] = min(leadDistance);
    closestVert(valsRight >= distanceThresh) = nan;
    leadToVertex(~isLeftInds) = closestVert;
    distanceFull(~isLeftInds) = valsRight;

end

leadLocationsPial = nan(numLeads,3);

isLeftInds = revertToVector(isLeftInds); 
assert(size(isLeftInds,1)==1 && size(leadToVertex,1)==1); 

leadLocationsPial(isLeftInds & ~isnan(leadToVertex),:) = surfL.vertices(leadToVertex(isLeftInds & ~isnan(leadToVertex)),:);
leadLocationsPial(~isLeftInds & ~isnan(leadToVertex),:) = surfR.vertices(leadToVertex(~isLeftInds & ~isnan(leadToVertex)),:);
% If use alternative is on, this maps STANDARD leadToVertex, using the ALT
% surface, to Euclidean space 

%% 

geodesicDistanceMaster = nan(numLeads,numVert);

parfor ii = 1:numLeads
    
    if isnan(leadToVertex(ii)); continue; end
    
    if isLeftInds(ii); geodesicDistanceMaster(ii,:) = myBp.dist_geodesic(surfL, leadToVertex(ii));
    else; geodesicDistanceMaster(ii,:) = myBp.dist_geodesic(surfR, leadToVertex(ii));
    end

end

geodesic.geodesicDistances = geodesicDistanceMaster;
geodesic.vertexNums = leadToVertex;
geodesic.isLeftInds = isLeftInds;
geodesic.leadLocationsPial = leadLocationsPial; 
geodesic.leadLocations = leadLocations;  % Dural 
geodesic.leadNames = leadNames; 
geodesic.distanceThresh = distanceThresh; 
geodesic.distanceFull = distanceFull; 
geodesic.useAlternative = useAlternative; 

fprintf('Geodesic distances built in %.3f seconds \n',toc)

if saving 
    if useAlternative; fn = fullfile('/Users/diamondjm/Documents/Josh_work/commonFiles/alternativeSurface/');
    else; fn = fullfile('/Users/diamondjm/Documents/Josh_work/commonFiles/');
    end
    if useUtah; fn = fullfile(fn,'geodesicUtah'); else; fn = fullfile(fn,'geodesic'); end
    fn = fullfile(fn,subj);
    save(fn,'geodesic');
end
    

end
