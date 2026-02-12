function badChans = retrieveBadChans(subj,sess)

badChansDir = '/Users/diamondjm/Documents/Josh_work/commonFiles/badChans/';


if exist(fullfile(badChansDir,subj,sess,'bad_chans.mat'),'file')
    badChans = load(fullfile(badChansDir,subj,sess,'bad_chans.mat'),'currentBadchans');
    badChans = badChans.currentBadchans;
else
    badChans = []; 
end


end

