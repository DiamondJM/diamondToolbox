function filenames = lsCell(directoryName, ISDIR)
% LSCELL returns the contents of ls as a cell array
%
% filenames = dirName(directoryName,[isdir])
%
% Equivalent to dir(directoryName).name --> filter out '.' names
% Does not return full paths, returns only the filenames
%
% Input
%   directoryName    - 
%   ISDIR (optional) - if given, filter results by returning only those satisfying isdir(filenames)==ISDIR
% 
    if nargin < 1 || isempty(directoryName)
        directoryName = pwd;
    end
    
    d = dir(directoryName);
    if nargin >= 2
        d = d([d.isdir] == ISDIR);
    end

    
    filenames = {d.name};
    
    mask = cellfun(@(x) x(1)~='.', filenames);
    filenames = filenames(mask);

end