function make_outer_surface (filled_volume, se_diameter, output_surface, dilation_radius)

% Original Author: Marie Schaer
% Date: 2007/11/14
% 
% New Author: Mike Trotta
% Date: 2018/07/11
%
% This function takes as an input the binary volume resulting of the
% filling of the any surface (usually the pial one) using mris_fill, and
% will close the sulci using morphological operation, with a sphere as the
% structural element.
% 
% Parameters: 
% se_diameter is the diameter of the sphere (in mm), use 10 mm by default 
% to close the sulci. After a point (usually 15), using a higher diameter has no effect
% 
% Utilities to write the surface to the freesurfer format are modified from
% "freesurfer_write_surf" (from Darren Weber's bioelectromagnetism toolbox), 
% according to a suggestion by Don Hagler in FreeSurfer's mailing list on 
% August 3, 2007.
%
%
% Example: make_outer_surface('lh.pial.mgz',15,'lh.outer-pial')
%
% REVISION HISTORY
%   07/18 MST - Logical code change was made to imdilate by 1 pixel after imclose
%             - Implementation change was made to keep image binary and use flat disks
%               in each of the x, y, and z directions instead of using ball as morphological
%               structure (for simplification + speedup)
%
% See Also: make_outer_surface

    if nargin < 4 || isempty(dilation_radius)
        dilation_radius = 2;
    end
    
    if isempty(se_diameter)
        se_diameter = 10;
    end
    
    diskc = strel('disk', se_diameter);
    diskd = strel('disk', dilation_radius);

    fprintf('reading filled volume...\n');
    temp = MRIread(filled_volume);
    im = logical(temp.vol);
    
    fprintf('closing volume...\n');
    imc = im;
    for i = 1:size(imc,1), imc(i,:,:)=imclose(squeeze(imc(i,:,:)), diskc); end
    for i = 1:size(imc,2), imc(:,i,:)=imclose(squeeze(imc(:,i,:)), diskc); end
    for i = 1:size(imc,3), imc(:,:,i)=imclose(squeeze(imc(:,:,i)), diskc); end
    
    fprintf('dilating volume...\n');
    imcd = imdilate(imc, diskd);
    
    fprintf('Creating outer surface...\n');
    [f,v] = isosurface(imcd);
    
    % 2D inspection
    DEBUG=0;
    if DEBUG
        dt = 0.05;
        figure();
        for i=1:size(im,1), imagesc(0.5*(squeeze(im(i,:,:)) + squeeze(imcd(i,:,:)))); title(i); pause(dt); end
        for i=1:size(im,2), imagesc(0.5*(squeeze(im(:,i,:)) + squeeze(imcd(:,i,:)))); title(i); pause(dt); end
        for i=1:size(im,3), imagesc(0.5*(squeeze(im(:,:,i)) + squeeze(imcd(:,:,i)))); title(i); pause(dt); end
    end
    
    % Mike did not make modifications below this point.
    
    v2=[129-v(:,1) v(:,3)-129 129-v(:,2)]; % in order to cope with the different orientation 
    v=v2;

    fprintf('Writing outer surface...\n');
    fname=output_surface;
    vert = v;
    face = f - 1;
    vnum = size(vert,1);  
    fnum = size(face,1);  
    
    % open it as a big-endian file
    fid = fopen(fname, 'wb', 'b');
    TRIANGLE_FILE_MAGIC_NUMBER = 16777214;
    fwrite3(fid, TRIANGLE_FILE_MAGIC_NUMBER);
    
    % Ouput a couple of text lines with creation date
    str = sprintf('created from matlab on %s\n',datestr(now));
    fwrite(fid, str,'char');
    fwrite(fid, vnum,'int32');
    fwrite(fid, fnum,'int32');

    % reshape vert into column array and write
    vert = reshape(vert',size(vert,1)*size(vert,2),1);
    fwrite(fid, vert,'float32');

    % reshape face into column array and write
    face = reshape(face',size(face,1)*size(face,2),1);
    fwrite(fid, face,'int32');
    fclose(fid) ;
end
