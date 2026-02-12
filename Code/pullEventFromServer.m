function [myStruct, failureFlag] = pullEventFromServer(subj,task)

[~,subj] = rectifySubj(subj); 

failureFlag = false; 
try
    myStruct = load(fullfile(pullServerDirectory, subj ,'/behavioral',task,'events.mat'));
    % myStruct = load(fullfile(pullServerDirectory(1), subj ,'/behavioral',task,'events.mat'));
    % warning('Old server directory.'); 
    myStruct = myStruct.events;
    
catch
    myStruct = struct; 
    failureFlag = true; 
    
end

end