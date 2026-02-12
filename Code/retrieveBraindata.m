function [myBd, myBp, useAlternative] = retrieveBraindata(subj,varargin)

%% Preamble 


p = inputParser;
addParameter(p,'forceNew',false);
addParameter(p,'saving',true); 
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
forceNew = p.Results.forceNew; 
saving = p.Results.saving; 

useAlternative = p.Results.useAlternative; 

% There was some weird stuff that went on with NIH023 and NIH046. Basically
% consider that there is the standard localized version (presumably Trotta
% made it). The issue is that the resection didn't work well with these
% patients. So Price re-made the localization. We'll call these the
% alternatives. They live in
% '/Users/diamondjm/Documents/Josh_work/commonFiles/AlternativeSurface/BD'.
% At one point I called them legacy, but that's confusing terminology. 
% They're legacy in that, if I remake all braindatas, the 'alternatives'
% get overwritten, and so I had to find them in a folder called legacy. 
% At any rate, Price's surfaces work well on the alternatives, but not on
% the standards. 
% It's hard to get a resection to fit the standard, so unfortunately the
% easier route is to simply localize using the alternative. This means that
% we need a version of geodesic that corresponds to the alternative
% anatomy. 
% UseAlternative can be passed to retrieveBraindata and loadGeodesic.

[~,subj] = rectifySubj(subj); % six digit

if ismember(subj,{'NIH023','NIH046'}) && ~useAlternative
    warning('Should use alternative be on for subject %s?',subj); 
    % useAlternative = true; 
end 


if useAlternative; baseDir = '/Users/diamondjm/Documents/Josh_work/commonFiles/AlternativeSurface/BD'; 
else; baseDir = '/Users/diamondjm/Documents/Josh_work/commonFiles/Braindata'; 
end

%% 

subjDir = fullfile(baseDir,subj); 

fName = fullfile(subjDir,'bdBpFull.mat');

if exist(fName','file') == 2 && ~forceNew 
    S = load(fName,'myBd','myBp');
    myBp = S.myBp;
    myBd = S.myBd;
else
    
    if forceNew && saving; rmdir(subjDir,'s'); end
    % This can be handled by determining whether or not forceNew is passed
    % into localizer_rois. 
    
    warning('Creating BDs de-novo. useFull is set to true. Note that, in order to get full lookup tables, we need to point to my local directory. So this change needs to be made in braindata''s constructors (lines 195-196)');
    % root = '/Volumes/Shares/FRNU/dataWorking/sz';
    root = pullServerDirectory; 
    myBd = braindata2(subj,root);
    myBp = ez_get_plotter(myBd);
    
    if saving
        mkdir(fullfile(subjDir));
        save(fullfile(subjDir,'bdBpFull'),'myBd','myBp')
    end
    
    
end

end