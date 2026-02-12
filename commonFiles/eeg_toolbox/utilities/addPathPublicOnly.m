function addPathPublicOnly(dirToAdd)
% Recursively adds all subdirectories with .m files
% to the path, excluding private directores [+@.'private']
% Depth-first recursion
%
% See also: MYPATH
    contents = dir(dirToAdd);
    hasMFolder = false;

    for i=1:length(contents)
        name = contents(i).name;
        contentPath = fullfile(dirToAdd,name);
        [~, ~, ext] = fileparts(contentPath);

        % recursively add non-private folders
        if (contents(i).isdir) && ...
                (~ismember(name(1),'.+@')) && ...   % exclude packages,methods,private dirs
                (~strcmp(name,'private'))           % exclude private dirs
            addPathPublicOnly(contentPath);   
            
        % mark this folder has .m files or class packages (@) or .mat files
        elseif strcmp(ext,'.m') || name(1) == '@' || strcmp(ext, '.mat')
            hasMFolder = true;
        end
    end

    if hasMFolder, addpath(dirToAdd), end;

end