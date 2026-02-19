function copyInterhemXYZ(filename, missingChanNames)
    % Designed to copy interhemispheric coordinates from one hemisphere to the other
    % Given a file with chanName, x, y, and z columns, and a cell array of new interhemispheric channels,
    % Add new rows for the missing channels with XYZ equal to the corresponding R/L-channels
    
    assert(exist(filename, 'file') > 0, 'Cant find %s', filename);
    t = readtable(filename);
    
    for i = 1:numel(missingChanNames)
        mischan = missingChanNames{i};
        
        if ismember(mischan, t.chanName)
            continue;
        end
        
        copyfromchan = mischan;
        if mischan(1) == 'L'
            copyfromchan(1) = 'R';
        elseif mischan(1) == 'R'
            copyfromchan(1) = 'L';
        else
            error('channel name %s expected to start with L or R', mischan);
        end
        
        assert(ismember(copyfromchan, t.chanName), 'Cant find an existing interhemispheric expected channel %s to match with %s', copyfromchan, mischan);
        
        irow = strcmpi(t.chanName, copyfromchan);
        row = t(irow, :);
        row.chanName{1} = mischan;
        t = [t; row];
        
        
    end
    
    writetable(t, filename);
    
end