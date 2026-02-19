function [res] = setBrainVizProps(res,hasColor)
if ~exist('hasColor','var')||isempty(hasColor)
  hasColor = false;
end

ax = gca;
ax.Clipping = 'off';

if ~hasColor
  set(res,'SpecularStrength',.1, ...
	  'DiffuseStrength',.6, ...
	  'SpecularColorReflectance',0, ...
	  'AmbientStrength',.45)
else
    set(res,'SpecularStrength',.1, ...
	  'DiffuseStrength',.4, ...
	  'SpecularColorReflectance',.5, ...
	  'AmbientStrength',.45)  
end
return