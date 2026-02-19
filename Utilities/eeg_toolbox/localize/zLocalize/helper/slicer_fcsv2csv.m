function anchor_table = slicer_fcsv2csv(filename)
    % Desc: Convert between fiducial csv and normal csv
    % Inputs:
    %   filename: path to the anchors.[f]csv file
    % Outputs:
    %   anchor_table

    [parent, name, ext] = fileparts(filename);
    filename_csv = [fullfile(parent, name) '.csv'];
    filename_fcsv = [fullfile(parent, name) '.fcsv'];
    
    if contains(filename, '.fcsv')
        conversion = 'TO_CSV';
    elseif contains(filename, '.csv')
        conversion = 'TO_FCSV';
    else
        error('file %s must be csv or fcsv', filename);
    end

    switch conversion
        case 'TO_CSV'
            anchor_table = slicer_read_fiducials(filename);
            anchor_table.useAnchor = ones(height(anchor_table),1);
            if exist(filename_csv, 'file')
                delete(filename_csv);
            end
            writetable(anchor_table, filename_csv);

        case 'TO_FCSV'
            anchor_table = readtable(filename);
            if ~ismember('chanName',anchor_table.Properties.VariableNames)
                anchor_table.chanName = util_combine_strnum(anchor_table.tagName, anchor_table.absElmtIdx);
                writetable(anchor_table, filename);
            end

            if exist(filename_fcsv, 'file')
                delete(filename_fcsv);
            end
            slicer_write_fiducials(anchor_table, filename_fcsv);
    end
end