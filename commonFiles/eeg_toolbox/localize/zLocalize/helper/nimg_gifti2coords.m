function my_coords = util_gifti2coords(gifti_path, varargin)
% Desc:
%   convert gifti to a coorrds table
% Inputs:
%   gifti_path: path to gifti file
% Outputs:
%   my_coords: coordinates for gifti 
%%
% read in gifti coordinates
my_gifti = gifti(gifti_path);
my_coords = double(my_gifti.vertices);

% outputs as xyz
my_coords = table(my_coords(:,1), my_coords(:,2), my_coords(:,3));
my_coords.Properties.VariableNames = {'x','y','z'};
return