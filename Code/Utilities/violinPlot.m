function violinPlot(x,y,d,color)

if isempty(y); return; end 

[k,xi] = ksdensity(y);
k = normalizeToBounds(k,[0 .4]); 
patch(k * d + x,xi,color)


end
