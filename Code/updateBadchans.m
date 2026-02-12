function currentBadchans = updateBadchans(subj,sess,newBadchans,varargin)

p = inputParser;
addParameter(p,'removeChans',false);
parse(p,varargin{:})
removeChans = p.Results.removeChans;

badChans = retrieveBadChans(subj,sess); 

if removeChans; currentBadchans = setdiff(badChans,newBadchans); 
else; currentBadchans = union(badChans,newBadchans); 
end

badChansDir = '/Users/diamondjm/Documents/Josh_work/commonFiles/badChans/';

dir = fullfile(badChansDir,subj); if ~isfolder(dir); mkdir(dir); end
dir = fullfile(dir,sess); if ~isfolder(dir); mkdir(dir); end
save(fullfile(badChansDir,subj,sess,'bad_chans.mat'),'currentBadchans');
        
end


