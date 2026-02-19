function n = count_filetype(d_parent, ext)
    % Counts the number of files in d_parent with the given extension
    % n = count_filetype(d_parent, ext)
    
    if ~iscellstr(ext), ext = cellstr(ext); end
    dir_struct  = dir(d_parent);
    names       = {dir_struct.name};
    [~,~,exts]  = cellfun(@fileparts, (names), 'uniformOutput',0); 
    mask        = ismember(upper(exts), strcat('.',upper(ext)));
    n           = sum(mask);
end