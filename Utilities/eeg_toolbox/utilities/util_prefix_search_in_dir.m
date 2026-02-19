function prefix_names = util_prefix_search_in_dir(my_dir, prefix, varargin)
% Desc:
%   find stuff with the allotted prefix in specified diretory
% Inputs:
%   my_dir
%   prefix
% Outputs:
%   prefix_names

%%
params.caseInsensitive = true; % dunno what happened here...
%%
% get directory contents
curdir = dir(my_dir);
dircon = {curdir.name};
dirlook = dircon;
fprintf('extracting files/directories...\n')
if params.caseInsensitive
    fprintf('\tcase insensitive...\n')
    dirlook = upper(dircon);
    prefix = upper(prefix);
end

% get the prefix names
prefix_idx= cellfun(@(x) ismember(1,x), strfind(dirlook, prefix),'UniformOutput',true);
prefix_names = dircon(prefix_idx);
fprintf('ahhhh sheeee...\n')
return