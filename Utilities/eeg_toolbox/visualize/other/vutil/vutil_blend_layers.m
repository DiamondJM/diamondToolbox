function blend_rgb = vutil_blend_layers(ulay_rgb, olay_rgb, blend)
% Desc:
% Inputs:
%   ulay_rgb: underlay_rgb
%   olay_rgb: overlay_rgb
%   blend: [0,1] strength of ulay 
% Outputs:
%   blend_rgb: blended rgb, assumes taht nans of olay are ignored 
% Notes: inspired by blendAnatomyData in figure_3d of brainstorm

% blend the shit that requires blending 
blend_rgb = ulay_rgb; % set overall as underlay, thus ignoring olay 
olay_vec = sum(olay_rgb,2);
to_blend = ~isnan(olay_vec);
blend_rgb(to_blend,:) = blend*ulay_rgb(to_blend,:) + (1-blend)*olay_rgb(to_blend,:);
return