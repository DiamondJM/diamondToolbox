% localizer_run_sample
% 
% This file runs the template localization pipeline.
% The directory structure will be: [rootDir]/ [subj]/ [subdirectories...]

rootDir = '~/zLocalize_data'; % <--- Edit this line - this is the folder where each subject will be stored
subj    = 'NIH046';           % this is the codename of the subject in sampleData
hemi    = 'r';                % subject had electrodes in a single hemisphere (right)

fprintf('\n\n---------------------------------------------------------------------------------------------------\n');
fprintf('Welcome to the zLocalize sample program. At this point you should have installed all necessary\n');
fprintf('dependencies (AFNI,Freesurfer,3DSlicer) and verified that they are in working condition, since\n');
fprintf('these programs will be called from the zLocalize MATLAB scripts. MATLAB can call other functions\n');
fprintf('in two different ways: (1) the folders containing those scripts (ie AFNI, Freesurfer, 3DSlicer)\n');
fprintf('can be on your system $PATH, and you can launch MATLAB directly from terminal so that MATLAB\n');
fprintf('inherits your system $PATH; or (2) You can edit the settings in localizer_run_template.m to specify\n');
fprintf('the full paths to the AFNI/Freesurfer/3DSlicer binaries\n');
fprintf('\n');
fprintf('You should also have added all subdirectories of zLocalize to your MATLAB path and downloaded the\n');
fprintf('sampleData folder. The sampleData folder should be separate from the data directories that zLocalize\n');
fprintf('will create. zLocalize stores all of its subects'' data in the rootDir directory. The first time it is\n');
fprintf('run for a subject, it will create its own directory structure of empty folders. It will then look for\n');
fprintf('all necessary input folders. If it does not find it in the expected location, it will launch a dialog\n');
fprintf('box which allows you to select the raw files it is looking for; it will copy those files into its\n');
fprintf('expected location\n');
fprintf('\n');

fprintf('** Read the READMEs in zLocalize and sampleData **\n');
fprintf('** When asked to locate files, navigate to your sampleData/ directory, and select the required files **\n');
fprintf('** For the "freesurfer subject directory," select sampleData/freesurfer/NIH046--not sampleData/freesurfer ** \n');
fprintf('** When asked for DICOM files, select Cancel **\n');
fprintf('** After your root data directory is created, you may copy over anchors.csv and anchors.fcsv from sampleData/ into\n');
fprintf('   [rootDir]/NIH046/anchors if you wish to forego placing your own anchors in Slicer (see README note) **\n');
fprintf('---------------------------------------------------------------------------------------------------\n');

input('Press [enter] to continue... ');
fprintf('\n\n');


if ~exist(rootDir,'dir'), mkdir(rootDir); end
localizer_run_template(subj, rootDir, hemi);