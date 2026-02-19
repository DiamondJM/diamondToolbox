function combined = surf_combine_hemis(both_surf)
combined = both_surf.lh; % initialize

for my_hemi = {'lh','rh'};
    both_surf.(my_hemi{:}).num_verts = length(both_surf.(my_hemi{:}).vertices);
end

combined.whichHemi = 'combined';
combined.vertices = [both_surf.lh.vertices; both_surf.rh.vertices];
combined.faces    = [both_surf.lh.faces; (both_surf.rh.faces + both_surf.lh.num_verts)];
combined.lhIdx = [1 both_surf.lh.num_verts];
combined.rhIdx = [(both_surf.lh.num_verts+1) (both_surf.lh.num_verts + both_surf.rh.num_verts)];
combined.path     = strrep(both_surf.lh.path, 'lh', '*');
combined.name     = strrep(both_surf.lh.path,'lh','*');
return