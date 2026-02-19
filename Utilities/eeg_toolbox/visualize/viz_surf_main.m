function [h_surf, surf_rgb] = viz_surf_main(my_surf, varargin)
    %   Below are a list of parameters
    %     'olay_rgb',  [], 'parameter';
    %     'ulay_rgb',  [], 'parameter';
    %     'olay_dset', [], 'parameter';
    %     'ulay_dset', [], 'parameter';
    %     'olay_cmap', struct(), 'parameter';
    %     'ulay_cmap', struct(), 'parameter';
    %     'olay_thresh',[NaN NaN],'parameter';
    %     'blend', 0, 'parameter'};
    %     'Tag'
    
    h_surf = [];
    % default ulay_cmap is a truncated gray scale map
    graymap = gray(100);
    graymap = graymap(90:-1:10,:);
    default_ulay_cmap.custom_cmap = graymap;

    default_cell = {'ulay_cmap',default_ulay_cmap,'parameter';
        'plotSurf',1, 'parameter'};
    params = util_extract_inputs(default_cell, varargin);
    
    num_pts = length(my_surf.vertices(:,1));
    params.num_pts = num_pts;

    surf_rgb = vutil_prep_layers(params);
    if params.plotSurf
        h_surf = vutil_view_surf(my_surf.vertices, my_surf.faces, surf_rgb, params);
    end
end

