function varargout = mypath()    
    % MYPATH lists the user's path without clutter of standard MATLAB paths
    % See also: STARTUP, path
    pathStr = path;
    pathCells = strsplit(pathStr,':');                 
    strfind(pathCells,matlabroot);                      
    isUserPath = cellfun(@isempty, strfind(pathCells, matlabroot));

    if nargout 
        varargout{1} = strjoin(pathCells(isUserPath),':');
    else
        disp(char(pathCells(isUserPath)));
    end
    
end


