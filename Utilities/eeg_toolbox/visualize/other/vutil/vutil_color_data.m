function [cdata cmap] = vutil_color_data(rdata, varargin)
% Desc: color my data!
% Inputs:
% Outputs
%%
if iscell(rdata)
    default_datatype = 'discrete';
else
    default_datatype = 'continuous';
end

%%
default_cell = {'datatype',default_datatype,'parameter';
    'name','parula','parameter'};
params = util_extract_inputs(default_cell, varargin);

%% color data based on data type
switch params.datatype
    case 'discrete'
        [cdata, cmap] = vutil_color_data_discrete(rdata, params);
    case 'continuous'
        [cdata, cmap] = vutil_color_data_continuous(rdata, params);
end
%%
return
function [cdata, cmap] = vutil_color_data_discrete(rdata, varargin)
% varargin = {}
%%
default_cell = {'custom_cmap',[],'parameter';
    'name','parula','parameter';
    'ordering','alphabetical','parameter'};
params = util_extract_inputs(default_cell, varargin);

%%
% sample data
% rdata = {'hi','bye','ay', 'hi'};

if isrow(rdata)
    rdata = rdata';
end

% get the uniques
uni_rdata = unique(rdata);
num_uni = length(uni_rdata);

% now create a map of that length 
cmap = table();
cmap.group = uni_rdata;
cmap.rgb = vutil_get_mlab_cmap(params.name, num_uni);

% uni 2 rgb table
map_table = table();
map_table.rdata = uni_rdata;
map_table.rgb = cmap.rgb;

rdata_table = table(rdata);

% join the tables to get cdata
rdata_table = join(rdata_table, map_table);
cdata = rdata_table.rgb;
%%

% %% 
% % get the color map that you want!
% if istable(cmap)
% else
%     if ~isempty(params.custom_cmap) && isstruct(params.custom_cmap)
%         cmap = custom_cmap;
%     elseif ~isempty(params.custom_cmap) && ~isstruct(params.custom_cmap)
%         cmap.clim = params.clim;
%         cmap.maps = params.maps;
%         cmap.rgb = custom_cmap;
%     else
%         cmap.clim = params.clim;
%         cmap.maps = params.maps;
%         cmap.rgb = vutil_get_mlab_cmap(params.name);
%     end
%     num_data = length(rdata);
%     cdata = cmapping(num_data, cmap, [1 1 1], 'colormap', 'discrete', cmaps.maps, cmap.clim);
% 
%     if strcmp(params.ordering,'alphabetical')
%         [sort_rdata, sidx] = sort(rdata);
%         cdata = cdata(sidx,:);
%     end
% end
%%
return
function [cdata, cmap] = vutil_color_data_continuous(rdata, varargin)
% Desc: 
%   for continuous variable. 
%   this thing will probably chagne someday, so don't pay too much
%   attention to it
% Inputs:
% Outputs:

% rdata = linspace(0,100);
% varargin = {'maps','direct'};
%%
default_cell = {'custom_cmap',[],'parameter';
    'name','parula','parameter';
    'clim',[nanmin(rdata) nanmax(rdata)],'parameter';
    'maps','scaled','parameter'};
params = util_extract_inputs(default_cell, varargin);

%%
if ~isempty(params.custom_cmap) && isstruct(params.custom_cmap)
    cmap = params.custom_cmap;
elseif ~isempty(params.custom_cmap) && ~isstruct(params.custom_cmap)
    cmap.clim = params.clim;
    cmap.maps = params.maps;
    cmap.rgb = params.custom_cmap;
else
    cmap.clim = params.clim;
    cmap.maps = params.maps;
    cmap.rgb = vutil_get_mlab_cmap(params.name, 64);
end
%% 
% sample data
nan_idx = isnan(rdata); % tabulate nans and put them back later

% carry out the mapping
cdata = cmapping(rdata, cmap.rgb, [1 1 1], 'colormap','continuous',cmap.maps, cmap.clim);
cdata(nan_idx,:) = nan(size(cdata(nan_idx,:)));
%%
% figure(2)
% scatter3(rdata, rdata, rdata, 30, cdata)
% axis off
%%
return
% function [cdata] = vutil_color_data(sdata, cmap)
% % Desc: 
% %   convert from scalar to color data
% % Inputs:
% %   |-dset: n x 1, data per vertex, roi, coordinate, etc.
% %   |-cmap: how to map from data to color
% %     |-Name: description YES?
% %     |-CMap: n_colors x 3 YES? 
% %     |-isAbsoluteValues: {0,1}: if 0 - scale, if 1 - absolute YES
% %     |-MaxValue: maximum value of color bar YES
% %     |-MinValue: minimum value of the colorbar YES
% %     |-Contrast: [-1, 1] - contrast for current colormap NOPE
% %     |-Brightness: [-1, 1] - Brightness for current colormap NOPE
% % Outputs:
% %   |-cdata: n x 3, rgb per vertex, roi, coordinate, etc. 
% 
% % if cmap is not provided, then provide one for us
% %   color map to sdata
% 
% return
