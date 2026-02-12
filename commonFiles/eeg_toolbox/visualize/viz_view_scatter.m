function [hgroup, hscat] = vutil_view_scatter(xyz, varargin)
% Desc:
%   Displays 3d scatter plot of target points
% Inputs:
%   xyz: num_pts x 3
%   varargin: 
%       dodge: 
%       size:
%       markerType:
%       filled: 
%       fillColor: 
% Outputs:
%   hgroup: handles group
%   hscat: collection of scatter handles
%% test case
% xyz = rand(100,3)
% % xyz = [0.5 0.5 0.5; 0 1 0; 1 0 0; 1 0.5 1];
% varargin = {};

%% set stuff up
num_pts = length(xyz(:,1));

%% define parameters
default_cell = {'dodge',[0 0 0], 'parameter';
    'size',1400,'parameter';
    'color',[1 0 0.5],'parameter';
    'markerType','o','parameter',;
    'tag','scatter','parameter';
    'filled',1,'parameter';
    'fillColor',[0.8 0.8 0.8],'parameter'};
params = util_extract_inputs(default_cell, varargin);

%% format parameters
expand_fields = {'dodge','size','color','markerType','tag','filled','fillColor'};
for my_name = expand_fields
    my_param = params.(my_name{1});
    if length(my_param(:,1)) < num_pts
        params.(my_name{1}) = repmat(my_param, num_pts, 1);
    end
end

xyz = xyz + params.dodge; % dodge the xyz
%% and display
% figure(); hold on;
for ipt = 1:num_pts
    hscat = scatter3(xyz(ipt,1), xyz(ipt,2), xyz(ipt,3), params.size(ipt), 'linewidth', 2, 'Marker', params.markerType(ipt),...
        'MarkerEdgeColor',params.color(ipt,:), 'MarkerFaceColor',params.fillColor(ipt,:));
%     hscat.inner(ipt) = scatter3(xyz(ipt,1),xyz(ipt,2),xyz(ipt,3),params.size(ipt), 'linewidth', 1, 'Marker',params.markerType(ipt));
%     if params.filled(ipt)
%         hscat.outer(ipt) = scatter3(xyz(ipt,1),xyz(ipt,2),xyz(ipt,3),params.size(ipt), params.fillColor(ipt,:), 'filled', 'linewidth',0.01);
%     else
%         hscat.outer(ipt) = [];
%     end
end
%%
% hgroup.inner = hggroup;
% hgroup.outer = hggroup;
hgroup = hggroup;

hgroup.Tag = params.tag(1);
set(hscat,'Parent', hgroup);
% set(hscat.outer,'Parent', hgroup.outer)
return