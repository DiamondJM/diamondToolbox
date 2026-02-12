%% ALICE pipeline
% % See https://github.com/UMCU-RIBS/ALICE
    % Local modifications by Mike here
    
% MP Branco et al. 2017
%
% (C) Mariana P. Branco, 
%     July 2017



%%%% Set paths:
close all; clear all; clc;

thisDirectory = pwd;

disp('                                         ');
disp('  ****** Welcome to ALICE ******          ');
disp('         Version July 2017 '                 );
disp('   A UMCU and NIH collaboration. ')
disp('                                         ');

addpathSystem('/Users/trottams/abin');


%%%% Some global GUI layout settings.
if ispc
	Settings.FontSize		= 8;
	Settings.GridFontSize	= 8;
elseif isunix
	Settings.FontSize		= 12;
	Settings.GridFontSize	= 10;
else
	Settings.FontSize		= 8;
	Settings.GridFontSize	= 8;
end

%Check dependencies:
%spmpath = which('spm');
system(['afni -help > afnipath.txt']);

f = fopen('afnipath.txt');
S = fscanf(f,'%s');

% if isempty(spmpath)
%     errordlg( 'Please add first SPM 12 to your path.','ERROR!');
% end
if isempty(S)
    errordlg( 'Please make sure to add AFNI and SUMA to your bash.','ERROR!');
    
else
    
    delete('afnipath.txt');
    
    %%%% Run program:
    gui = ctmrGUI;
    
end




