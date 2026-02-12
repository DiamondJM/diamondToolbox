function startTime = readEEG_wrapper(subj,sess,task)

% This should work only for NK. 


rawPath = fullfile(pullServerDirectory,subj,'raw',task,sess); 
rawFiles = lsCell(rawPath); 

rawPrefix = cell(size(rawFiles)); 
for ii = 1:length(rawFiles)
    [~,rawPrefix{ii}] = fileparts(rawFiles{ii});
    if contains(rawPrefix{ii},'21E')
        [~,rawPrefix{ii}] = fileparts(rawPrefix{ii});
    end
end

assert(length(unique(rawPrefix)) == 1); 

eegFile = fullfile(pullServerDirectory,subj,'raw',task,sess,sprintf('%s.EEG',rawPrefix{1})); 

[timestamp,T_second] = readEEG(eegFile); 

startTime = datetime(timestamp,'InputFormat','yyMMdd_HHmm') + seconds(T_second);

end