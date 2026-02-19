function view_table = vutil_format_coords(coords_table, varargin)
% Desc:
% Inputs:
% Outputs:
%   view_table: ripe and ready for viewing
%%
default_cell = {
    'color','','parameter';
    'colorOpts',struct(),'parameter';
    'fill','','parameter';
    'fillOpts',struct(),'parameter';
    'text','','parameter';
    'size','','parameter';
    'sizeOpts',struct(),'parameter'
    };
inputs = util_extract_inputs(default_cell,varargin);
%% prepare stuff for viewing
view_table = table();
if ~ismember('xyz',coords_table.Properties.VariableNames)
    view_table.xyz = coords_table{:,{'x','y','z'}};
else
    view_table.xyz = coords_table.xyz;
end

% both color and fill opts have the options of
% custom_cmap
% name
% clim
% maps
% ordering
% datatype

if ~isempty(inputs.color)
    color_me = coords_table.(inputs.color);
    view_table.color = vutil_color_data(color_me, inputs.colorOpts);
end

if ~isempty(inputs.fill)
    fill_me = coords_table.(inputs.fill);
    view_table.fill = vutil_color_data(fill_me, inputs.fillOpts);
end

if ~isempty(inputs.text)
    view_table.text = coords_table.(inputs.text);
end

if ~isempty(inputs.size)
    view_table.size = coords_table.(inputs.size);
end
return