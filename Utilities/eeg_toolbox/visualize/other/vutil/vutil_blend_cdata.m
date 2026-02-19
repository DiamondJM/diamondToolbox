function [mixrgb] = vutil_blend_cdata(urgb, orgb, toBlend, ualpha)
% Input:
%   |-urgb:   overlay rgb
%   |-orgb:   underlay rgb
%   |-toBlend: blend these vertices
%   |-ualpha: strength with which to view ulay
% Outputs:
%   |-mixrgb: mixed color data
% Credit:
%   lifted code from 'BlendAnatomyData' in brainstorm
mixrgb = urgb; 
if ~isempty(orgb)
    mixrgb(toBlend,:) = ualpha * urgb(toBlend,:) + (1 - ualpha) * orgb(toBlend,:);
end
return
