function serverDir = pullServerDirectory(useOld)

if nargin == 0; useOld = false; end

if ~useOld; serverDir = '/Volumes/joshSzExtraction/eeg'; 
else; serverDir = '/Volumes/Josh_myBook/szBackup/sz';
end

% serverDir = '/Volumes/Shares/FRNU/dataWorking/sz'; 

end