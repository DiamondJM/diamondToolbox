function all_grow_table = mesh_grow_centers(my_surf, center_list, radius, log)
    % Desc: wrapper to grow geodesic circles about centers
    % Inputs:
    %   log - name of file to log to (append). If not passed, logged to stdout
    % Outputs:
    %   all_grow_table with columns:
    %       ROIC_idx      - Index 1 to k of this row's ROI center among all k ROI centers
    %       ROIC_mesh_idx - Index into surface mesh of the ROIC (in [1,n])
    %       vertex        - Index into pial surface mesh (standard SUMA-resampled)
    %       d             - Distance from ROIC to vertex (mm)
    %       
    %       

    
    % previous table column names were changed from the following:
    %   grow_struct = struct('to_list',nan,'d',nan,'vert_idx',nan,'center_idx',nan);
    
    grow_struct = struct('vertex',nan,'d',nan,'ROIC_mesh_ndx',nan,'ROIC_ndx',nan);
    grow_struct_ar = repmat(grow_struct, length(center_list), 1);

    n2f = surfing_nodeidxs2faceidxs_custom(my_surf.faces',my_surf.vertices');
    V = double(my_surf.vertices);
    F = my_surf.faces;
    face_verts = my_surf.faces(:);

    %parfor i_center = 1:length(center_list)
    for i_center = 1:length(center_list)
        ROIC_mesh_ndx = center_list(i_center);

        missing = ~ismember(ROIC_mesh_ndx, face_verts);
        if ~missing        

            [vertex , d, ~, ~] = surfing_circleROI_custom(V, F, ROIC_mesh_ndx, radius, 'geodesic',n2f);

        grow_struct_ar(i_center).vertex = vertex';
        grow_struct_ar(i_center).d = d';
        grow_struct_ar(i_center).ROIC_mesh_ndx = repmat(ROIC_mesh_ndx,length(vertex),1);
        grow_struct_ar(i_center).ROIC_ndx = repmat(i_center,length(vertex),1);


        end
    end

    % Consolidate
    vertex = cat(1, grow_struct_ar.vertex);
    d = cat(1, grow_struct_ar.d);
    ROIC_mesh_ndx = cat(1, grow_struct_ar.ROIC_mesh_ndx);
    ROIC_ndx = cat(1, grow_struct_ar.ROIC_ndx);
    all_grow_table = table(ROIC_ndx, ROIC_mesh_ndx, vertex, d);

end
