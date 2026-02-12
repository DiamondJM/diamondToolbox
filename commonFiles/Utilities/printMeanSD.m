function outputString = printMeanSD(valList)

meanVal = mean(valList,'omitnan');
stdVal = std(valList,'omitnan');

% warning('Using confidence interval.'); 
% stdVal = calculateConfInterval(valList);
outputString = sprintf('%.2f %s %.2f',meanVal,'\pm',stdVal'); 

fprintf('%s\n',outputString); 

end