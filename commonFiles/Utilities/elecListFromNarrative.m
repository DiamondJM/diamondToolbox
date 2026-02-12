function outList = elecListFromNarrative(subj)

[~,subj] = rectifySubj(subj); 
masterList = generateChansMasterList(subj); 

outList = {}; 

while true 

    tn = input('Please enter tag name. Type q to quit. ','s'); 
    if isequal(tn,'q'); break; end

    n = input('Please enter numbers. '); 

    for ii = 1:length(n)
        outList{end + 1} = sprintf('%s%d',tn,n(ii)); 
    end

    assert(all(ismember(outList,masterList))); 
end

outList = unique(outList); % Also sorts

fn = fullfile('leadsMappedPositive',sprintf('%s.mat',subj)); 
save(fn,'outList'); 

end
