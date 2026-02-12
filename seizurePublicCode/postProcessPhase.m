function [dataMaster,params] = postProcessPhase(dataMaster)


%%%%%%
% Now post-process this data in such a way as to only pick the good data
%%%%%

% quantileVal = 0.10;

phase = dataMaster.phase;
pow = dataMaster.pow; 


quantileVal = 0; 
quantileThresh = quantile(pow(:),quantileVal);

retainChans = 14;

for jj = 1:size(phase,1)
    
    currentPow = pow(jj,:);
    [~, ind] = sort(currentPow,'descend');
    
    lowPowInds = ind(retainChans + 1:end);
    
    phase(jj,lowPowInds) = nan;
end

phase(pow < quantileThresh) = nan;

pow(isnan(phase)) = nan; 

% cycleLength = 1 / freqsRange(1) * Fs / stepSize;
% cycleLength = 1 / min(filterWindow) * Fs / stepSize;


cycleLength = nan; 
% cycleLength = 100;
% for ii = 1:size(phaseMaster,1)
%     x = [0 ~isnan(phaseMaster(ii,:)) 0];
%     idx1 = strfind(x,[1 0]) - 1;
%     idx0 = strfind(x,[0 1]);
%     [sortLength, sortInds] = sort(idx1 - idx0 + 1,'descend');
%     
%     sortInds = sortInds(sortLength > cycleLength);
%     
%     preserveLog = false(size(phaseMaster(ii,:)));
%     
%     for iii = sortInds
%         preserveLog(idx0(iii) : idx1(iii)) = true;
%     end
%     
%     phaseMaster(ii, ~preserveLog) = nan;
% end


fprintf('%.f percent of all phase data used. \n',sum(~isnan(phase(:))) / numel(phase) * 100);

dataMaster.phase = phase;
dataMaster.pow = pow;


params.cycleLength = cycleLength; 
params.quantileVal = quantileVal; 
params.retainChans = retainChans;