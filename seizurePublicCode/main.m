%% Localization pipeline. 

% Hit run to use our pipeline for seizure source localization with an
% example structure corresponding to a single seizure for a single patient. 

addpath(genpath(pwd));
load('myStruct.mat','myStruct')

%% Source localization

myStruct = localizationManager(myStruct); 

%% Plotting

figure; 
plotSurfFun(myStruct,0); 

figure; 
plotDimRedVert1(myStruct); 

figure; 
plotDimRedVert2(myStruct); 