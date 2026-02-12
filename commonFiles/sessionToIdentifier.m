function identifier = sessionToIdentifier(mySingleton)

subj = mySingleton.subjID; 

currentDate = mySingleton.date;
currentDate = datestr(currentDate,'yymmdd-HHMMss'); 

identifier = sprintf('%s_%s',subj,currentDate); 

end