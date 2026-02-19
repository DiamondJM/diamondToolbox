function [blend_rgb, olay_cmap, ulay_cmap] = vutil_prep_layers(varargin)
    % Desc: prepare the surface rgb 
    % Inputs:
    %   see default_cell...
    % Outputs:
    %   blend_rgb: the blend_rgb to use with the shiz
    %% extract parameters/inputs
    default_cell = {
        'num_pts',   [], 'parameter';
        'olay_rgb',  [], 'parameter';
        'ulay_rgb',  [], 'parameter';
        'olay_dset', [], 'parameter';
        'ulay_dset', [], 'parameter';
        'olay_cmap', struct(), 'parameter';
        'ulay_cmap', struct(), 'parameter';
        'olay_thresh',[NaN NaN],'parameter';
        'blend', 0, 'parameter'};
    params = util_extract_inputs(default_cell, varargin);


    %% when either overlay or underlay are not provided, output SOMETHING
    if isempty(params.ulay_rgb) && isempty(params.ulay_dset) % color the underlay if nothing provided
        params.ulay_rgb = repmat([0.5 0.5 0.5], params.num_pts, 1);
    end

    if isempty(params.olay_rgb) && isempty(params.olay_dset) % turn off overlay with NaNs if nothing provided
        params.olay_rgb = repmat([NaN NaN NaN], params.num_pts, 1);
    end

    %% and now handle all of the information
    % prepare underlay and overlay rgb
    % Mike's notes: the below functions only resort to rgb if there's no dset. Really there isn't a good way to mix rgb and data
    %               Either the surface has its color property set to a value+cmap or to an rgb
    [olay_rgb, olay_cmap] = vutil_process_layer('rgb', params.olay_rgb, 'dset',params.olay_dset, 'thresh', params.olay_thresh, 'cmap' ,params.olay_cmap, 'num_pts', params.num_pts);
    [ulay_rgb, ulay_cmap] = vutil_process_layer('rgb', params.ulay_rgb, 'dset', params.ulay_dset, 'cmap' ,params.ulay_cmap,'num_pts',params.num_pts);

    %% blend it all together
    % blend the overlay and underlay, where blend = 0 means overlay takes full
    % color priority over a vertex
    blend_rgb = vutil_blend_layers(ulay_rgb, olay_rgb, params.blend);
    %%
end