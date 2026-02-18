function plotAllDirections(myStruct)

lLat = [270,0];
rLat = [-270,0];
vInf = [90,-90];
vSup = [90,90];

clf
subplot(2,2,1)
myStruct.plotSurfFun;
view(lLat(1),lLat(2))
ax = gca;

copyobj(ax.Children,subplot(2,2,2))
myBp.apply_axis_settings();
myBp.apply_surf_settings();
view(rLat(1),rLat(2))

copyobj(ax.Children,subplot(2,2,3))
myBp.apply_axis_settings();
myBp.apply_surf_settings();
view(vInf(1),vInf(2))

copyobj(ax.Children,subplot(2,2,4))
myBp.apply_axis_settings();
myBp.apply_surf_settings();
view(vSup(1),vSup(2))


% [~,subj] = rectifySubj(patients(ii));
% saveas(gcf,fullfile('/Volumes/USB DISK/Josh/',subj,sprintf('sz_%d.svg',kk)))


end