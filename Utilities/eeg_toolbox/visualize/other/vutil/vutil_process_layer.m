function [rgb cmap] = vutil_process_layer(varargin)
% Desc:
%   take in either a dset + cmap or an rgb dset
%   doesn't really do anything if it is an rgb dset
%   if dset, then colors the dset... so, really, it's just 
%   vutil_color_data
% Inputs:
%   dset + cmap
%   rgb
% Outputs:
%   rgb
%   cmap
% Notes:
%   dunno if this funciton is really necessary

% set defaults
default_cell = {'rgb',[],'parameter';
    'dset',[],'parameter';
    'cmap',struct(),'parameter'; % inherits default cmap
    'thresh',[NaN NaN],'parameter';
    'num_pts',1,'parameter'}; 
inputs = util_extract_inputs(default_cell, varargin);

% switch modes based on what is provided
if ~isempty(inputs.rgb) && ~isempty(inputs.dset)
    error('only input rgb or dset + cmap!');
elseif ~isempty(inputs.dset)
    dset = inputs.dset;
    
    % handle thresholding if provided
    if ~isnan(inputs.thresh(1))
        dset(dset < inputs.thresh(1)) = NaN;
    end
    if ~isnan(inputs.thresh(2))
        dset(dset > inputs.thresh(2)) = NaN;
    end
    
    % convert dset to rgb
    rgb = vutil_color_data(dset, inputs.cmap);
    cmap = inputs.cmap;
elseif ~isempty(inputs.rgb)
    % prepare if only one color provided
    if length(inputs.rgb(:,1)) == 1
        rgb = repmat(inputs.rgb, inputs.num_pts, 1);
    else
        rgb = inputs.rgb;
    end
    cmap = [];
end
return 