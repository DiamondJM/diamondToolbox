function [hgroup, hsphere] = vutil_view_sphere(xyz, varargin)
    % Desc:
    %   Displays 3d spheres
    % Inputs:
    %   xyz: num_pts x 3
    %   varargin: name-value pair
    %       radius: 
    %       color:
    %       dodge:
    %       tag:
    %       figure:
    % Outputs:
    %   hgroup: handles group
    %   hsphere: collection of surf handles

    %% test case
    % xyz = [0 0 1; 0 1 0; 1 0 0; 1 1 1];
    % varargin = {};

    %% set stuff up
    num_pts = length(xyz(:,1));
    N_SPHERE = 10;

    ip = inputParser;
    ip.addParameter('radius', 1.5);
    ip.addParameter('color', [1 1 1]);
    ip.addParameter('tag', 'electrode');
    ip.addParameter('opacity', 1);
    ip.addParameter('figure', []); % (Cocjin's, Mike doesn't use this)
    ip.addParameter('dodge', [0 0 0]); % (Cocjin's, Mike doesn't use this)
    ip.KeepUnmatched = 1;
    ip.parse(varargin{:});
    

    

    % Replicate parameters
%     results = ip.Results;
%     fns = fieldnames(ip.Results);
%     for fni = numel(fns)
%         param = ip.Results.(fns{fni});
%         if length(param(:,1)) < num_pts
%             results.(fns{fni}) = repmat(param, num_pts, 1);
%         end
%     end
    results = ip.Results;
    
    color = results.color;
    
    s = table();
    [s.x, s.y, s.z] = sphere(N_SPHERE); % sphere object

   
    xyz = xyz + ones(size(xyz,1),1)*ip.Results.dodge; % dodge the xyz  %JW added matrix multiplication to make dimensions match (perhaps just required for older matlab version

    
    for ipt = 1:num_pts
        sph = s;
        i_dodge = 1;
        for my_ax = {'x','y','z'}
            sph.(my_ax{1}) = results.radius * s.(my_ax{1}) + xyz(ipt,i_dodge);
            i_dodge = i_dodge + 1;
        end

        if size(color, 2) == 3
            % RGB
            if size(color,1) == 1
                c = color;
            else
                c = color(ipt,:);
            end

            hsphere(ipt) = surf(sph.x, sph.y, sph.z, ...
                            'FaceColor',c,...
                            'EdgeColor','none');
        else
            % data (colormap)
            if numel(color) == 1
                c = repmat(color, N_SPHERE+1, N_SPHERE+1);
            else
                c = repmat(color(ipt), N_SPHERE+1, N_SPHERE+1);
            end

            hsphere(ipt) = surf(sph.x, sph.y, sph.z, c, ...
                            'FaceColor','flat',...
                            'EdgeColor','none');
        end

        hold on;
        set(hsphere(ipt),'FaceLighting','phong',...
            'AmbientStrength',0.5,...
            'SpecularStrength', 0, ...
            'FaceAlpha',results.opacity,...
            'Tag',results.tag);
        % create a new surface object. the surface data is the scaled unit
        % sphere where the center has been translated by the coordinate value
        % of the data point that the sphere represents.
    end
    %%
    hgroup = hggroup;
    set(hsphere,'Parent',hgroup);
    hgroup.Tag = results.tag;
    %%
end