%% lifted from Episurg
%% Set Lighting & View
shading interp; lighting gouraud; material dull; axis off, hold on
if ischar(brainView)
    switch brainView
        case 'r'
            l=light('Position',[1 0 0]);
            view(90,0)
        case 'rm'
            l=light('Position',[-1 0 0]);
            view(270,0)
        case 'rim'
            l=light('Position',[-1 0 0]);
            view(270,-45)
        case 'ri'
            l=light('Position',[0 0 -1]);
            view(90,-90)
        case 'ro'
            l=light('Position',[0 -1 0]);
            view(0,0)
        case 'lo'
            l=light('Position',[0 -1 0]);
            view(0,0)
        case 'rf'
            l=light('Position',[0 1 0]);
            view(180,0)
        case 'lf'
            l=light('Position',[0 1 0]);
            view(180,0)
        case 'rs'
            l=light('Position',[0 0 1]);
            view(90,90);
        case 'rsv' %superior & vertically aligned
            l=light('Position',[0 0 1]);
            view(0,90);
        case 'l'
            l=light('Position',[-1 0 0]);
            view(270,0);
        case 'lm'
            l=light('Position',[1 0 0]);
            view(90,0);
        case 'li'
            l=light('Position',[0 0 -1]);
            view(90,-90);
        case 'lim'
            l=light('Position',[-1 0 0]);
            view(270,-45);
        case 'ls'
            l=light('Position',[0 0 1]);
            view(270,90);
        case 'lsv' %superior & vertically aligned
            l=light('Position',[0 0 1]);
            view(0,90);
        case 'liv' %inferior & vertically aligned
            l=light('Position',[0 0 -1]);
            view(0,-90);
        case 'riv' %inferior & vertically aligned
            l=light('Position',[0 0 -1]);
            view(0,-90);
    end
    clear l
else
    light('Position',brainView.light);
    view(brainView.eyes)
end
alpha(opaqueness);