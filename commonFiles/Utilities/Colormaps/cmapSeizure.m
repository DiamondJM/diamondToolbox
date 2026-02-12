function cmap = cmapSeizure

cmap = [interp1([0;1],[[68 57 107] / 256; [159 149 199] / 256],linspace(0,1,128));interp1([0;1],[[159 149 199] / 256; [221 212 252] / 256],linspace(0,1,128))];
% cmap = flipud(cmap); 

end