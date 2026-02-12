function [new_face_list] = surf_trim_faces(face_list, vert_subset, varargin)
% Desc:
%   weed them faces
% Inputs:
%   face_list: n_face x 3
%   vert_subset: subset of them faces
% Outputs:
%   new_face_list: the surviving faces

%% parse variable inputs
ip = inputParser;
ip.addParameter('overlap_num',2);
ip.parse(varargin{:});

%% weed out them faces and return them 
face_idx = sum(ismember(face_list,vert_subset),2) > ip.Results.overlap_num;
new_face_list = face_list(face_idx,:);
return
