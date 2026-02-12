function dataVals = normalizeToBounds(dataVals,endRange,beginRange) 


if nargin <= 2 
    beginRange = [min(dataVals) max(dataVals)]; 
end
if nargin <= 1 
    endRange = [0 1];
end

%% Put between 1 and 0 

dataVals = max(beginRange(1),min(dataVals,beginRange(2)));
% Normalize to bounds
if range(beginRange) > 0
    dataVals = (dataVals - beginRange(1)) / range(beginRange);
else
    dataVals = zeros(size(dataVals)); 
end

%% Put into end range 

dataVals = dataVals * range(endRange) + endRange(1); 



