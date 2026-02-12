% visualization_utahs.m 
% 
% Author: Mike Trotta - adapted from Cocjin's visualization_script.m


clear; clc; close all;


%%
subjs_utahs = [29,30,34,37,39,42];
subj = 'NIH034';
rootEEGdir = '/Volumes/Shares/FRNU/dataWorking/eeg/FINAL_RUN';
surfType = 'pial';

path_focus = 'best';
params_cell = {'subjId' subj;
               'nodeSurf','std141';
               'parcType','zroi'};
path_list = util_flood_list(params_cell, util_list_paths('rootEEGdir', rootEEGdir, 'path_focus', path_focus));
dir_list = util_flood_list(params_cell, util_list_dirs('rootEEGdir', rootEEGdir, 'path_focus', path_focus));
jacksheet = readtable(path_list.jack);
utahs = jacksheet(strcmpi('UTAH', jacksheet.hardwareType),:);
coords = readtableSafe(fullfile(dir_list.monopolar, 'compiled.csv'));
utahs = join(utahs,coords);
xyz = utahs{:,{'x','y','z'}};
% note: these coords^ are on the individual's dural surface

%%
% now we have coordinate and surfaces
% snap these coordinates to the standard pial surface
bd_pia = braindata('subjId', subj, 'suma_dir', dir_list.suma, 'dural', false);
bd_pia.fn_load_surf('surfType','pial');

%% Find the nearest std.141 vertex on this subjects's pia
[nearestVertex, dist] = knnsearch(bd_pia.surf.vertices, xyz);
utahs.nearestVertex = nearestVertex;
%bd_pia = [];

%%
color = [0 1 0]; % green
utahs.color = repmat(color, height(utahs), 1);

%%
% fg = viz_build_grid(1, [1 3]);
% fg(1) = figure();
clear bd
hemispheres = {'lh','rh'};
for i = 1 : length(hemispheres), whichHemi = hemispheres{i};
    % load the subject's surface up
    bd.(whichHemi) = braindata('subjId',subj,'suma_dir', dir_list.suma);
    bd.(whichHemi).fn_load_surf('surfType',surfType, 'whichHemi', whichHemi);
    
    figure('name', bd.(whichHemi).meta.subjId);
    curv = bd.(whichHemi).curv.data;
    curv(curv < 0) = -1;
    curv(curv > 0) = 1 ;
    a = viz_surf_main(bd.(whichHemi).surf, 'ulay_dset',curv);
    viz_set_shine(a);
    axes(a.Parent); %#ok<LAXES>
    axis equal; shading interp; lighting gouraud; material dull; axis off; hold on;
    view(235,0)
    
    
    %
    utahs_hem = strcmp(utahs.whichHemi, whichHemi);
    surf_idx = utahs{utahs_hem, 'nearestVertex'};
    %
	xyz = bd.(whichHemi).surf.vertices(surf_idx,:);
    if ~isempty(xyz)
        viz_view_sphere(xyz, 'color', utahs{utahs_hem,'color'}, 'fill', utahs{utahs_hem,'color'}, 'radius',2 )
    end
    %
end
camlight;

