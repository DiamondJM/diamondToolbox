function patientInfo = parsePatientInfo(subj, rootEEGdir, resection_string_only)
%
%  Load a patient info file and parse out the various bits of info contained within
%
%
% REVISION HISTORY
%   04/2018 MST - use getJackTable and add resection_string_only parameter

if nargin < 3
    resection_string_only = 0;
end



patientInfoFile = fullfileEEG(rootEEGdir,subj,'docs/patient_info.txt');

subjNum = str2num(subj(4:6)); %- NIHXXX or BEHXXX


%- prep the return structure (even if it cant open the file)
clear patientInfo;
patientInfo.subject       = ''; %-placeholder
patientInfo.subjNum       = subjNum;
patientInfo.filePath      = patientInfoFile;
patientInfo.ANY_TAG_MISSING = 1;
patientInfo.TAG_COUNT_GOOD  = 0;
patientInfo.TAG_COUNT_TOTAL = 0;




%- confirm patient info file exists and can be opened
if ~exist(patientInfoFile, 'file')
    error('Cant find file: %s', patientInfoFile); 
end

if resection_string_only
    try
        output = '';
        cmd = sprintf('grep -i resection %s', patientInfoFile);
        [~,output] = system(cmd);
        patientInfo.resection = strtrim(output);
        return
    catch e
        fprintf('Error searching for resection string: %s\n', output);
        fprintf('Error: %s\n', e.message);
    end
end

fid = fopen(patientInfoFile,'r');
if fid == -1
    error('Problem opening file: %s', patientInfoFile); 
end

%[subj, rootEEGdir] = fileparts(subjectPath);
jacktable = getJackTable(subj, rootEEGdir);
chanTags = jacktable.chanName;
[~, chanNums] = util_split_stringnum(chanTags);


%- add in the jacksheet info
patientInfo.jacksheetTags = chanTags;
patientInfo.jacksheetNums = chanNums;



%-  Example data:
% Subject: NIH0XX
% Age at test: 33
% Ethnicity: white
% Gender: male
% Handedness: right
% language dominance: left
% resection:


searchParams   = {'Subject:', 'Age at test:','Ethnicity:','Gender:','Handedness:','language dominance:','resection:'};  %- this is the expected string in the patient_info file
saveParams     = {'subject',  'ageAtTest',   'ethnicity', 'gender', 'handedness', 'langaugeDominance',  'resection'};   %- this will be the field names in the output struct


breakString    = '------------------';
searchClinical = {'Ictal:', 'Interictal:', 'Resected:', 'IsCut:'};  %- all text between these strings and the break string are assumed to be electrodes
saveClinical   = {'Ictal',  'Interictal',  'Resected',  'IsCut'};  %- all text between these strings and the break string are assumed to be electrodes



%- initialize variables
fileText     = '';
fileTextCell = {};
paramVals    = cell(size(saveParams));
foundParams  = zeros(1,length(searchParams));
foundClinic  = zeros(1,length(searchClinical));


%- next read the file and parse the "searchParams" lines.  save file contents to cell array for further parsing below
while 1, 
    tline = fgetl(fid);
    
    %- break when you find the end of the important info
    %if length(tline)==0 | ~ischar(tline), break; end;
    if ~ischar(tline), break; end;
    
    %- save a copy of the entire file
    if length(fileText)==0, fileText=tline; else fileText = sprintf('%s\n%s',fileText,tline); end
    fileTextCell{end+1} = tline;
    
    
    %- also parse out specific info
    for iParam=1:length(searchParams),
        thisString = searchParams{iParam};
        if length(tline)>=length(thisString) && strcmp(thisString,tline(1:length(thisString))),
            paramVals{iParam} = strtrim(tline(length(thisString)+1:end)); %- string trim cuts leading and training spaces
            foundClinic(iParam) = 1;
        end
    end
end
if ~contains(fileTextCell{end},breakString(1:15)), fileTextCell{end+1} = breakString; end %- and a terminating breakstring if not there already... make code below more robust


%- close the file when done
fclose(fid);


%- populated the basic output structure
patientInfo.fileText     = fileText;
patientInfo.fileTextCell = fileTextCell;
for iParam = 1:length(searchParams),
    patientInfo.(saveParams{iParam}) = paramVals{iParam};
end



%-------------------------------------------------------------------------------
%---- populate electrode-based clinical info

%- first create list of all important indicies so that beginning and end of each category can be identified
%iAll = find(strcmp(fileTextCell,breakString));  %- specific number of characters here... need to enforce this
iAll = find(contains(fileTextCell,breakString(1:15)));  %- relax requirement for exact number of dashes... now just 15 dashes or more equals a break.

%if length(iAll)==0, fprintf('\n heads up, no hyphen dividers found... looking for exactly %s',breakString); end
for iParam = 1:length(searchClinical)
    iAll = [iAll find(strcmp(fileTextCell,searchClinical{iParam}))]; %- complete list of relevant indicies into fileTextCell
    %patientInfo.(sprintf('%s_grouped',saveClinical{iParam})) = {};  %- placeholder
end
iAll = unique([iAll length(fileTextCell)]); %- include end of file if not there already


ANY_TAG_MISSING = 0;
TAG_COUNT_GOOD  = 0;
TAG_COUNT_TOTAL = 0;


%- identify electrodes within each category
for iParam = 1:length(searchClinical)
    
    thisParam  = searchClinical{iParam};
    paramStart = find(contains(fileTextCell,thisParam));
    
    %- assume not listed... starting around NIH044 these are only present for "none"
    iList = 0;  chanList = {};  isNone=0;
    
    if isempty(paramStart), 
        if ~strcmp(thisParam,'IsCut:') & (subjNum<=42 & strcmp(subj(1:3),'NIH')),
            %- isCut is optional, so only give the warning for the other ones
            fprintf('\nUh oh.. didnt find header %s in NIH subject before 43',thisParam);
        end
        %chanList = {'nan 0'}; %- do we want a dummy channel when nothing is written down?
    else
        %- cell array with entries between "Ictal" and "Interictal".... could be blank or contain "None"
        potentialChan = fileTextCell(paramStart+1:iAll(min(find(iAll>paramStart)))-1);
        
        %- create a cell array just with channel lists (e.g. RALT 1-4)... still need to convert to list of individual channels
        while iList<length(potentialChan),
            iList=iList+1;
            tline = potentialChan{iList};
            if     length(tline)==0, continue;
            elseif length(tline)>=3 && strcmp(tline(1:3),'---'), warning('\n WARNING: looks like hyphen-break exists, but is different length than expected'); keyboard;
            elseif length(tline)>=4 && strcmp(lower(tline),'none'), chanList = {}; isNone=1; break;
            else   chanList{end+1} = tline;
            end
        end
    end
    
    %- save a binary that conveys if patient info says "none"... if so, then we dont expect to find any electrodes
    patientInfo.(sprintf('%s_isNone',saveClinical{iParam})) = isNone;
    
    %- save copy of version direct from text (grouped by strip)
    patientInfo.(sprintf('%s_grouped',saveClinical{iParam})) = chanList;
    
    
    %- now split each group into a list of individual channels, and confirm they map to jacksheet master
    chanSplit   = {};
    chanJackNum = [];
    for iList=1:length(chanList),
        tline = chanList{iList};
        
        iSpace   = find(isspace(tline),1,'first');
        thisTag  = tline(1:iSpace-1);
        
        
        %- parse the channel list
        C = textscan(tline(iSpace+1:end),'%d%c');
        iNum = 1;   numList = [];
        while iNum<=length(C{1}),
            if length(C{2})<iNum, %- last number in a list (would already be caught) or single number in a list (only caught here)
                numList = [numList C{1}(iNum)];
                iNum = iNum+1;
            elseif C{2}(iNum)=='-',
                numList = [numList C{1}(iNum):C{1}(iNum+1)];
                iNum = iNum+1; %- skip the next number
            elseif C{2}(iNum)==','
                numList = [numList C{1}(iNum) C{1}(iNum+1)];
                iNum = iNum+1;
            else
                error('\n shouldnt be here');
            end
        end
        numList = unique(numList); %- oversample the list above then trim it here
        %disp(numList);
        
        
        %- now convert to tagnames and confirm a match up with jacksheet master
        for iNum=1:length(numList),
            thisChanTag = sprintf('%s%d',thisTag,numList(iNum));
            chanSplit{end+1} = thisChanTag;
            chanJackNum(end+1) = nan;
            
            iJacksheet = find(strcmp(chanTags,thisChanTag));
            if     length(iJacksheet)==0, fprintf('\n uh oh. constructed channel (%s) doesnt exist in jacksheet', thisChanTag); ANY_TAG_MISSING=1; %disp(chanTags);       %- this could happen if patient info had a typo OR if channel "isNotRec" but was actually there
            elseif length(iJacksheet)>1,  fprintf('\n uh oh. constructed channel (%s) has >1 entry in jacksheet. saving first match', thisChanTag); ANY_TAG_MISSING=1;   %- this should never happen
            else   chanJackNum(end) = chanNums(iJacksheet);  TAG_COUNT_GOOD = TAG_COUNT_GOOD+1;
            end
        end
    end
    
      
    %- save copy of split versions (each channel separate... this is the real goal)
    patientInfo.(sprintf('%s_split',saveClinical{iParam}))   = chanSplit;
    patientInfo.(sprintf('%s_chanNum',saveClinical{iParam})) = chanJackNum;
  
    TAG_COUNT_TOTAL = TAG_COUNT_TOTAL+length(chanSplit);
end


patientInfo.ANY_TAG_MISSING = ANY_TAG_MISSING; 
patientInfo.TAG_COUNT_GOOD  = TAG_COUNT_GOOD;
patientInfo.TAG_COUNT_TOTAL = TAG_COUNT_TOTAL;





