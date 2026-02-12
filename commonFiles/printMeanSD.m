function outputString = printMeanSD(valList)

meanVal = mean(valList,'omitnan');
stdVal = std(valList,'omitnan');

% warning('Using confidence interval.'); 
% stdVal = calculateConfInterval(valList);

% outputString = sprintf('%.1e %s %.1e',meanVal,'\pm',stdVal'); 
outputString = sprintf('%.2f %s %.2f \n',meanVal,'pm',stdVal); 
fprintf(outputString)

end