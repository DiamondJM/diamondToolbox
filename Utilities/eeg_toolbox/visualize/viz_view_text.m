function [hgroup, htext] = vutil_view_text(xyz, my_text, varargin)
% Desc:
%   Displays text in 3d space
% Inputs:
%   xyz: num_pts x 3
%   my_text: num_pts x (text or num)
%   varargin: name-value pair
%       dodge: 
%       size: 
%       color: 
%       tag: 
% Outputs:
%   hgroup: handles group
%   hgroup: collection of text handles
%% test case
% xyz = [0 0 1; 0 1 0; 1 0 0; 1 1 1];
% my_text = [0 1 2 3];
% varargin = {};

%% set up stuff
num_pts = length(xyz(:,1));
if isnumeric(my_text)
    my_text = arrayfun(@num2str, my_text, 'UniformOutput',false);
end

%% define parameters
default_cell = {'dodge', [0 0 0] ,'parameter';
    'size', 14, 'parameter';
    'color', [0.5 0.5 0.5],'parameter';
    'tag','text','parameter'
    };
params = util_extract_inputs(default_cell, varargin);

%% format parameters
for my_name = fieldnames(params)'
    my_param = params.(my_name{1});
    if length(my_param(:,1)) < num_pts
        params.(my_name{1}) = repmat(my_param, num_pts, 1);
    end
end
xyz = xyz + params.dodge;
%% and display
for ipt= 1:num_pts
    htext(ipt) = text(xyz(ipt,1), xyz(ipt,2), xyz(ipt,3), my_text(ipt));
    set(htext(ipt),...
        'FontSize',params.size(ipt,:),...
        'Color',params.color(ipt,:), ...
        'FontWeight','Bold',...
        'HorizontalAlignment','center',...
        'Tag',params.tag(ipt,:));

end
hgroup = hggroup;
set(htext,'Parent',hgroup);
hgroup.Tag = params.tag(1);
%%
return