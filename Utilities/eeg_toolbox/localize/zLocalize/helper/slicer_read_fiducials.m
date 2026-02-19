function leadsInfo = slicer_read_fiducials(fid_path, varargin)
    % copy into proper csv file
    [pathstr,name,~] = fileparts(fid_path);
    fcsv_path = fullfile(pathstr,[name '.fcsv']);
    csv_path  = fullfile(pathstr,[name '.csv']);
    copyfile(fcsv_path, csv_path);

    % read in as table
    fid_table = readtable(csv_path,'HeaderLines',3,'ReadVariableNames',0);
    delete(csv_path);

    var_name_list = {'id','x','y','z','ow','ox','oy','oz','vis','sel','lock','label','desc','associatedNodeID'}; % by default
    if isempty(fid_table)
       fid_table = cell2table(cell(0,length(var_name_list))); 
    end
    fid_table.Properties.VariableNames = var_name_list;

    % convert to leads table
    leadsInfo = fid_table(:,{'label','x','y','z'});
    leadsInfo.Properties.VariableNames = {'chanName','x','y','z'};

    % extract out tagName, absElmtIdx, and tagMod
    chanName = leadsInfo{:,'chanName'};
    isDig = arrayfun(@(x) isstrprop(x{:},'digit'),chanName,'UniformOutput',false);
    tagName = {};
    absElmtIdx = [];
    tagMod = {};
    for iLead = 1:length(chanName)
        % read in for lead
        iTagIdx = chanName{iLead};
        numChar = length(iTagIdx);
        iDigLog = find((isDig{iLead}));

        % extract tagName and absElmtIdx
        if ~isempty(iDigLog)
            tagName{iLead,1} = iTagIdx(1:(iDigLog(1)-1)); %#ok<*AGROW>
            absElmtIdx(iLead,1) = str2double(iTagIdx(iDigLog));
            tagMod{iLead,1} = iTagIdx(iDigLog(end)+1:numChar);
        else
            tagName{iLead,1} = iTagIdx;
            absElmtIdx(iLead,1) = NaN;
            tagMod{iLead,1} = '';
        end


    end
    leadsInfo.tagName = tagName;
    leadsInfo.absElmtIdx = absElmtIdx;
    leadsInfo.tagMod = tagMod;
    leadsInfo.orientation = repmat({'LPI'},height(leadsInfo),1);
end
