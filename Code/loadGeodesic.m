function geodesic = loadGeodesic(sourceLoc,varargin)

%% Preamble

p = inputParser;
addParameter(p,'unModified',false);
addParameter(p,'forceNew',false);
addParameter(p,'useAlternative',false); 
addParameter(p,'useUtah',false);
addParameter(p,'saving',true);
parse(p,varargin{:})
unModified = p.Results.unModified;
forceNew = p.Results.forceNew;
useAlternative = p.Results.useAlternative;
useUtah = p.Results.useUtah;
saving = p.Results.saving; 

subj = sourceLoc.subj;
[~,subj] = rectifySubj(subj);

if ismember(subj,{'NIH023','NIH046'}) && ~useAlternative 
    warning('Should use alternative be on for subject %s?',subj); 
    % useAlternative = true; 
end 

%% Load or create, unpack 

if useAlternative; fn = fullfile('/Users/diamondjm/Documents/Josh_work/commonFiles/alternativeSurface/');  
else; fn = fullfile('/Users/diamondjm/Documents/Josh_work/commonFiles/'); 
end
if useUtah; fn = fullfile(fn,'geodesicUtah'); else; fn = fullfile(fn,'geodesic'); end 
fn = fullfile(fn,sprintf('%s.mat',subj)); 
    
if exist(fn,'file') == 2 && ~forceNew
    load(fn,'geodesic')
else
    geodesic = buildGeodesicDistances_master(subj,'useUtah',useUtah,'useAlternative',useAlternative,'saving',saving);
end


%% Modify

if unModified; return; end

%% Unpack


chanNames = sourceLoc.chanNames;
[chanLog,chanInds] = ismember(chanNames,geodesic.leadNames);
% missingLeads = setdiff(chanNames,leadsFromGeodesic);
missingLeads = chanNames(~chanLog); 

if ~all(chanLog)
    warning('Channels submitted that are missing from our records.');
    
    disp(missingLeads); 
    % JD 5/6/2023: 
    % The new (improved) paradigm is that geodesicDistances should be
    % created ONLY using leads.csv. In other words, functionality which
    % cross-checks leads with time series (typically Jacksheet leads) with
    % leads with localization, has been removed. We should use only and all
    % the leads in leads.csv to create geodesic. Then pull relevant leads
    % from that 'master' list later. 
    
    % So as of that time, there's no solid reason to programatically
    % re-make geodesic. 
        
end



% The above doesn't work anymore because we need to permit tolerance of
% missing values. 
% Under this new paradigm (3/6/23) we will proceed with localization even
% if there are channels in jacksheet that failed to localize. 

leadsFromGeodesic = repmat({''},size(chanNames)); 
leadsFromGeodesic(chanLog) = geodesic.leadNames(chanInds(chanLog)); 
% Should match chanNames, minus the missing channels 

geodesicDistances = nan(length(chanNames),size(geodesic.geodesicDistances,2)); 
geodesicDistances(chanLog,:) = geodesic.geodesicDistances(chanInds(chanLog),:); 

vertexNums = nan(1,length(chanNames)); 
vertexNums(chanLog) = geodesic.vertexNums(chanInds(chanLog)); 

isLeftInds = false(size(chanNames)); 
isLeftInds(chanLog) = geodesic.isLeftInds(chanInds(chanLog)); 
% A bit risque to initialize this is as false. It should be fine, though,
% because all other fields are nan so we won't be able to do anything with
% the unlocalized ones (chanNames(~chanLog)). 

[leadLocations,leadLocationsPial] = deal(nan(length(chanInds),3)); 
leadLocations(chanLog,:) = geodesic.leadLocations(chanInds(chanLog),:); 
leadLocationsPial(chanLog,:) = geodesic.leadLocationsPial(chanInds(chanLog),:); 

distanceFull = nan(1,length(chanNames));
distanceFull(chanLog) = geodesic.distanceFull(chanInds(chanLog));

%% Pack

geodesic.geodesicDistances = geodesicDistances;
geodesic.vertexNums = vertexNums;
geodesic.isLeftInds = isLeftInds;
geodesic.leadLocations = leadLocations;
geodesic.leadLocationsPial = leadLocationsPial;
geodesic.leadNames = leadsFromGeodesic;

geodesic.distanceFull = distanceFull;

end
