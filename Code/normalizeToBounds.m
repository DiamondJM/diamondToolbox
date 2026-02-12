function dataVals = normalizeToBounds(dataVals,endRange,beginRange) 


wasFlipped = false; 
if size(dataVals,2) > 1 && size(dataVals,1) == 1 % Row vector
    dataVals = dataVals'; 
    wasFlipped = true; 
end

for ii = 1:size(dataVals,2)

    dataCurrent = dataVals(:,ii);


    if nargin <= 2 || all(isnan(beginRange))
        beginRange = [min(dataCurrent) max(dataCurrent)];
    end
    if nargin <= 1
        endRange = [0 1];
    end

    %% Put between 1 and 0

    dataCurrent = max(beginRange(1),min(dataCurrent,beginRange(2),'includenan'),'includenan');
    % Normalize to bounds
    if range(beginRange) > 0
        dataCurrent = (dataCurrent - beginRange(1)) / range(beginRange);
    else
        dataCurrent = dataCurrent - beginRange(1);
    end

    %% Put into end range

    dataCurrent = dataCurrent * range(endRange) + endRange(1);

    dataVals(:,ii) = dataCurrent;
end

if wasFlipped; dataVals = dataVals'; end

