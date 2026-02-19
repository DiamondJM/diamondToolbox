function [hgroup, hcyl] = viz_view_cylinder(xyz, pial_surf, varargin)
% [hgroup, hcyl] = vutil_view_cylinder(xyz, varargin) Plots cylinder(s) on the current axis
% Mike made this function and based it on Cocjin's viz_view_sphere
%   
% The cylinder is oriented according to the average surface normal in a small area around
% the given xyz-coordinate
%
% Input:
%   xyz     - xyz of the center of 1 face of the cylinder (e.g. electrode xyz coordinate); n X 3 matrix
%   pial_surf- struct with vertices and faces fields (e.g. a single pial hemisphere surface)
%   
% Optional Input:
%   radius  - radius of cylinder (default 1.5)
%   depth   - depth of cylinder (default 1)
%   color   - single color or an n X 3 RGB (1 row for each cylinder)
%   tag     - tag of graphics object to label cylinders (default 'electrode')
%       
% Output:
%   hgroup  - handle to group of cylinder objects
%   hcyl    - handle to clinder objects (graphics object array)
%
% Revision History
%   02/2018 MST - Created
%
% To Do:
%   - take into account grid/strip neighbors (include in a weighted surface normal average)
%   - load normal-region from file instead of on-the-fly calculation 

    N_CYL = 16;         % lower for redering speed, higher for quality (drawCylinder default is 32)
    R_CONTACT = 1.5;    % surface-normal region will be at least this big
    R_GROW = 1;          % surface-normal region will be padded with this
    
    
    
    
    ip = inputParser;
    ip.addParameter('radius', 1.5);
    ip.addParameter('color', [0 1 0]);
    ip.addParameter('tag', 'electrode');
    ip.addParameter('depth', 1);
    ip.addParameter('figure', []); % (Cocjin's, Mike doesn't use this)
    ip.addParameter('dodge', [0 0 0]); % (Cocjin's, Mike doesn't use this)
    ip.KeepUnmatched = 1;
    ip.parse(varargin{:});
    results = ip.Results;    
    
    color = results.color;
    depth = results.depth;
    num_pts = length(xyz(:,1));
    
    F = double(pial_surf.faces);
    V = double(pial_surf.vertices);
    n2f = surfing_nodeidxs2faceidxs(F');
   
    % Find the small area around the given coordinate
    lead_mesh_lut = lead_to_mesh_grow(pial_surf, R_CONTACT, R_GROW, xyz);
    
    for ipt = 1:num_pts


        % turn it into an actual subsurface
        sv = lead_mesh_lut.geodisc_mesh_ind{ipt};
        [sv, sf] = surfing_subsurface(V, F, sv, [], n2f);
        %p = patch('faces',sf, 'vertices',sv);
        

        % Mean of surface normals ( in the correct direction)
        norms = patchnormals(sv,sf);
        pnorm = mean(norms);
        pnorm = pnorm ./ norm(pnorm);
        xyz0 = xyz(ipt,:);
        if dot(pnorm, xyz0) < 0
            pnorm = -pnorm;
        end

        % create the cylinder
        xyz1 = xyz0 + pnorm * depth;
        cyl = cat(2, xyz0(:)', xyz1(:)', R_CONTACT);
        
        
        % draw
        if size(color, 2) == 3
            % RGB
            if size(color,1) == 1
                c = color;
            else
                c = color(ipt,:);
            end

            hcyl(ipt) = drawCylinder(cyl, N_CYL, 'closed', 'FaceColor', color);
        else
            % data (colormap)
            if numel(color) == 1
                c = repmat(color, N_CYL+1, N_CYL+1);
            else
                c = repmat(color(ipt), N_CYL+1, N_CYL+1);
            end

            hcyl(ipt) = drawCylinder(cyl, N_CYL, 'closed', 'FaceColor', 'flat');
        end
        
        set(hcyl(ipt),  'Tag', results.tag)
                        


    end
    
    hgroup = hggroup;
    set(hcyl,'Parent',hgroup);
    hgroup.Tag = results.tag;
    
end

function N = patchnormals(sv, sf) 
    %Vertex normals of a triangulated mesh, area weighted, left-hand-rule 
    % N = patchnormals(FV) -struct with fields, faces Nx3 and vertices Mx3 
    %N: vertex normals as Mx3

    %face corners index 
    A = sf(:,1); 
    B = sf(:,2); 
    C = sf(:,3);

    %face normals 
    n = cross(sv(A,:)-sv(B,:),sv(C,:)-sv(A,:)); %area weighted

    %vertice normals 
    N = zeros(size(sv)); %init vertix normals 
    for i = 1:size(sf,1) %step through faces (a vertex can be reference any number of times) 
        N(A(i),:) = N(A(i),:)+n(i,:); %sum face normals 
        N(B(i),:) = N(B(i),:)+n(i,:); 
        N(C(i),:) = N(C(i),:)+n(i,:); 
    end
end