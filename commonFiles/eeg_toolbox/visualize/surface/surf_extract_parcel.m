function sub_surf = surf_extract_parcel(my_surf, labels, vert_list, targ, varargin)
    % extract parcel... or the negative of the parcel

    %% parse input variables
    ip = inputParser;
    ip.addParameter('inverse',0);
    ip.parse(varargin{:});

    %% this is probably not the best way to do this
    vert_idx = surf_query_verts(labels, vert_list, targ);
    num_verts = length(my_surf.vertices(:,1));

    % go from vert_idx to sub_vert_idx to handle the inverse case as well as 
    % for repeat vert_idx - kinda rough
    zlog = zeros(num_verts, 1);
    zlog(vert_idx) = 1;
    if ip.Results.inverse
        zlog = ~zlog;
    end
    sub_vert_idx = find(zlog);

    sub_surf = my_surf;
    sub_surf.parc_vert_idx = sub_vert_idx;
    sub_surf.faces = surf_trim_faces(my_surf.faces, sub_vert_idx);
    %%
end

