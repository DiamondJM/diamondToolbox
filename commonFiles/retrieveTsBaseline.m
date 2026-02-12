function tsBaseline = retrieveTsBaseline(mySingleton,lengthBaseline)

error('Deprecated.');
% use populateTs. 

if nargin == 1; lengthBaseline = 30; end

currentDate = generateDate(mySingleton);
currentDate = sprintf('ts_%s.mat',currentDate);
 
filename = fullfile(sprintf('/Users/diamondjm/Documents/Josh_work/commonFiles/tsBaseline%d/',lengthBaseline),currentDate);
if exist(filename,'file') == 2 
    load(filename,'tsBaseline')

else
    currentSubj = mySingleton.subjID;
    currentEpoch = mySingleton.eegfile(end-10:end);
    startTime = mySingleton.eegoffset / 1000 - lengthBaseline;
    endTime = lengthBaseline;
    chanNames = mySingleton.sourceLoc.chanNames;
    
    tsBaseline = ...
        load_eegfile(currentSubj, currentEpoch, startTime,endTime,...
        'loc_detrend', 1, 'rmline', 3, 'rem_sat', 10, ...
        'filter', 'lowpassJosh', 'visual',0,'hardware','subdural','channels',chanNames);

    % save(filename,'tsBaseline')
    
end
