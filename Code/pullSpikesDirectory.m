
function spikesDirectory = pullSpikesDirectory(useOld)

if nargin == 0; useOld = false; end 

if ~useOld; spikesDirectory = '/Users/diamondjm/Documents/Josh_work/SourceLocSpikes/2023Feb25';
else; spikesDirectory = '/Users/diamondjm/Documents/Josh_work/SourceLocSpikes/2022Apr10'; warning('Old.'); 
end
% spikesDirectory = '/Users/diamondjm/Documents/Josh_work/SourceLocSpikes/2022Feb6';
% spikesDirectory = '/Users/diamondjm/Documents/Josh_work/SourceLocSpikes/2022Jan26';

end