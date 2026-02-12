function viewMRI(subj)

[~,subj] = rectifySubj(subj); 

% currentFolder = pwd; 

addpath(genpath('/Users/diamondjm/Documents/MATLAB/Utilities/fieldtrip-20180123'));
ft_defaults

image = ft_read_mri(fullfile(pullServerDirectory,subj,'/tal/zloc/freesurfer',subj,'mri/brain.mgz'));
% image = ft_read_mri(fullfile(sprintf('/Volumes/Shares/FRNU/data/eeg/%s/tal/zloc/CT_1/ct_implant.nii',subj)));

% image = ft_read_mri(fullfile(pullServerDirectory,subj,'/tal/zloc/mr_pre','mr_pre.nii')); 

clf; % axis image; 
colormap bone 
axis equal 
for ii = length(image.anatomy):-1:1
    imagesc(rot90(squeeze(image.anatomy(:,ii,:)),3));
    % if ii == 1; 
    set(gca,'ydir','normal');
    axis equal
    pause(.05);
end

% cd(currentFolder);

% preOvernightChecklist;

end