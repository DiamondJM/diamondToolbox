function cmap = vutil_get_mlab_cmap(cmap_name, n)
switch cmap_name
    case 'parula'
        cmap = parula(n);
    case 'jet'
        cmap = jet(n);
    case 'winter'
        cmap = winter(n);
    case 'cool'
        cmap = cool(n);
    case 'spring'
        cmap = spring(n);
    case 'summer'
        cmap = summer(n);
    case 'autumn'
        cmap = autumn(n);
    case 'hsv'
        cmap = hsv(n);
    case 'hot'
        cmap = hot(n);
    case 'gray'
        cmap = gray(n);
    case 'bone'
        cmap = bone(n);
    case 'copper'
        cmap = copper(n);
    case 'pink'
        cmap = pink(n);
    case 'rb'
    otherwise
        error('I do not recognize cmap of type: %s',cmap_name);
end
return