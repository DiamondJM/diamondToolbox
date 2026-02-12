function filepath = uigetfile_if_not_exist(default_path, filespec, title, sym, isdir, isBatch)

    if nargin < 3, title = ''; end
    if nargin < 4, sym = 0; end
    if nargin < 5, isdir = 0; end
    if nargin < 6, isBatch = 0; end
    
    if exist(default_path, 'file') && ~(isdir && is_empty_dir(default_path))
        filepath = default_path;
    else
        if isBatch
            filepath = [];
        else
            fprintf('%s\n', title);
            if isdir
                filepath = uigetdir(default_path, title);
            else
                [filename,pathname] = uigetfile(filespec, title, fileparts(default_path));
                filepath = fullfile(pathname, filename);
            end
        end
        if ischar(filepath) && exist(filepath, 'file') && ~strcmpi(filepath, default_path)
            if sym
                if exist(default_path, 'dir')
                    rmdir(filepath);
                end
                fprintf('Linking %s <-- %s\n', filepath, default_path);
                unix(sprintf('ln -s %s %s', filepath, default_path));
            else
                fprintf('Copying %s --> %s\n', filepath, default_path);
                copyfile(filepath, default_path);
            end
        else
            filepath = [];
        end
    end
end

function empty = is_empty_dir(filepath)
    try
        count = numel(lsCell(filepath));
    catch
        s = dir(filepath);
        charname = char(s.name);
        count = sum(charname(:,1) ~= '.');
    end
    empty = count == 0;
end