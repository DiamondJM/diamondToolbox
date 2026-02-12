function [masterList,hasSubdural,hasDepth,leadsSubdural,leadsDepth] = generateChansMasterList(subj)
% function [masterList,hasSubdural,hasDepth,leadsSubdural,leadsDepth] = generateChansMasterList(subj,rootEEGDir)

% if nargin < 2; rootEEGDir = pullServerDirectory; end

rootEEGDir = pullServerDirectory;

[~,subj] = rectifySubj(subj); 

leadsSubdural = getLeads(subj,rootEEGDir,'hardwareType','subdural'); 
leadsDepth = getLeads(subj,rootEEGDir,'hardwareType','depth'); 

hasSubdural = ~isempty(leadsSubdural); hasDepth = ~isempty(leadsDepth); 
 
masterList = [leadsSubdural; leadsDepth]; 

end