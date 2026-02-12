classdef ctmrGUI < handle
    % See https://github.com/UMCU-RIBS/ALICE
    % Local modifications by Mike here
    
    properties
        mainFig
        extraFig
        controls
        settings
        subGUI
    end
    
    
    methods
        
        function obj = ctmrGUI
            % Create main window.
            width		= 1000;
            height		= 750;
            frameheight = 200;
            
            %number lines in action log.
            obj.settings.NUM_LINES = 9;
            
            %no files loaded.
            obj.settings.loaded = [0 0 0];
            
            %initial values for 3dclustering
            obj.settings.CV = 3999;
            obj.settings.R  = 3;
            obj.settings.IS = 1;
            
            %default hemisphere:
            obj.settings.Hemisphere = 'Left';
            
            %default method:
            obj.settings.Method = 'Trotta';
            
            %empty grid settings:
            obj.settings.Grids = [];
            
            % Get screen size.
            screenSize = get(0,'ScreenSize');
            
            % Starting point of each frame.
            startingPointFrame = [10 round(width/3)+10 2*round(width/3)+10];
            
            % Main window
            windowPosition = [ round((screenSize(3)-width)/5), screenSize(4)-height-100, width, height+80];
            obj.mainFig = figure( 'Name', 'ALICE: Also Developed by a Cast of Thousands','OuterPosition', windowPosition, 'Menu', 'none', ...
                'NumberTitle', 'off', 'Color', get(0,'DefaultUIControlBackgroundColor'), 'Resize', 'off', 'CloseRequestFcn', @obj.figCloseRequest );
            
            %two buttons for create directory or locate directory:
            obj.controls.btnCreateDirectory = uicontrol( 'Parent', obj.mainFig, 'Style', 'pushbutton', 'Position', [295 718+40 200 35], ...
                'String', 'Create Directory', 'Callback', @obj.btnCreateDirectory, 'FontSize', 13 , 'FontWeight', 'bold');
            
            obj.controls.btnLocateDirectory = uicontrol( 'Parent', obj.mainFig, 'Style', 'pushbutton', 'Position', [505 718+40 200 35], ...
                'String', 'Locate Directory', 'Callback', @obj.btnLocateDirectory, 'FontSize', 13 , 'FontWeight', 'bold');
            
            %Frame 1
            obj.controls.frame1 = uipanel( 'Parent', obj.mainFig, 'Units', 'pixels', 'Position', [startingPointFrame(1) frameheight+40 round(width/3)-20 height-250], ...
                'Title', '1. CT-MR Co-registration', 'FontSize', 14,'Visible', 'off', 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [.5 .5 .5] );
            
            %Frame 2
            obj.controls.frame2 = uipanel( 'Parent', obj.mainFig, 'Units', 'pixels', 'Position', [startingPointFrame(2) frameheight+40 round(width/3)-20 height-250], ...
                'Title', '2. Electrode Selection', 'FontSize', 14, 'Visible', 'off','FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [.5 .5 .5] );
            
            %Frame 3
            obj.controls.frame3 = uipanel( 'Parent', obj.mainFig, 'Units', 'pixels', 'Position', [startingPointFrame(3) frameheight+40 round(width/3)-20 height-250], ...
                'Title', '3. Electrode Projection', 'FontSize', 14,'Visible', 'off', 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [.5 .5 .5] );
            
            %% Inside frame 1:
            % SubFrame 1 inside frame 1
            obj.controls.subframe1 = uipanel( 'Parent', obj.controls.frame1, 'Units', 'pixels', 'Position', [10 360 293 100], ...
                'Title', 'Select MRI scan', 'FontSize', 12, 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [0.8 .8 .8] );
            
            % SubFrame 2 inside frame 1
            obj.controls.subframe2 = uipanel( 'Parent', obj.controls.frame1, 'Units', 'pixels', 'Position', [10 240 293 100], ...
                'Title', 'Select FreeSurfer segmentation', 'FontSize', 12, 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [0.8 .8 .8] );
            
            % SubFrame 3 inside frame 1
            obj.controls.subframe3 = uipanel( 'Parent', obj.controls.frame1, 'Units', 'pixels', 'Position', [10 120 293 100], ...
                'Title', 'Select CT scan', 'FontSize', 12, 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [0.8 .8 .8] );
            
            % Button 1 inside frame 1
            obj.controls.btnAlignCTtoMRI = uicontrol( 'Parent', obj.controls.frame1, 'Style', 'pushbutton', 'Position', [90 35 140 50], ...
                'String', 'Align CT to MRI', 'Callback', @obj.btnAlignCTtoMRI, 'FontSize', 11, 'FontWeight', 'bold' );
            
            %% Inside frame 2:
            % SubFrame 1 inside frame 2
            obj.controls.subframe4 = uipanel( 'Parent', obj.controls.frame2, 'Units', 'pixels', 'Position', [10 360 293 100], ...
                'Title', 'CT scan selected in Step 1', 'FontSize', 12, 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [0.8 .8 .8] );
            
            % SubFrame 2 inside frame 2
            obj.controls.subframe5 = uipanel( 'Parent', obj.controls.frame2, 'Units', 'pixels', 'Position', [10 120 293 220], ...
                'Title', '3D-Clustering settings', 'FontSize', 12, 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [0.8 .8 .8] );
            
            % Button 1 inside frame 2
            obj.controls.btnExtractClusters = uicontrol( 'Parent', obj.controls.frame2, 'Style', 'pushbutton', 'Position', [10 35 140 50], ...
                'String', 'Extract Clusters', 'Callback', @obj.btnExtractClusters, 'FontSize', 11 , 'FontWeight', 'bold');
            
            % Button 2 inside frame 2
            obj.controls.btnSelectElectrodes = uicontrol( 'Parent', obj.controls.frame2, 'Style', 'pushbutton', 'Position', [163 35 140 50], ...
                'String', 'Select Electrodes', 'Callback', @obj.btnSelectElectrodes, 'FontSize', 11 , 'FontWeight', 'bold');
            
            
            %% Inside frame 3 ========================================================================================================================
            % SubFrame 1 inside frame 3
            obj.controls.subframe6 = uipanel( 'Parent', obj.controls.frame3, 'Units', 'pixels', 'Position', [10 120 293 340], ...
                'Title', 'Select Projection Method', 'FontSize', 12, 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [0.8 .8 .8] );
            
            % Button 1 inside frame 3
            obj.controls.btnVisualize = uicontrol( 'Parent', obj.controls.frame3, 'Style', 'pushbutton', 'Position', [90 35 140 50], ...
                'String', 'Visualize!', 'Callback', @obj.btnVisualize, 'FontSize', 11 , 'FontWeight', 'bold', 'Enable', 'off');
            
            
            %% Inside Subframe 1 - frame 1  =================================================
            % Button 1 inside subframe 1
            obj.controls.btnOpenMRI = uicontrol( 'Parent', obj.controls.subframe1, 'Style', 'pushbutton', 'Position', [10 25 100 40], ...
                'String', 'Open', 'Callback', @obj.btnOpenMRI, 'FontSize', 11 , 'FontWeight', 'bold');
            
            % text box inside subframe 1
            obj.controls.txtMRI = uicontrol( 'Parent', obj.controls.subframe1, 'Style', 'edit', 'Position', [115 26 168 36], ...
                'FontSize', 12, 'string', {' MRI scan (*.nii)'} , 'HorizontalAlignment', 'left', 'BackgroundColor', 'w' ,'enable','inactive');
            
            %% Inside subframe 2 - frame 1 =================================================
            % Button 1 inside subframe 2
            obj.controls.btnOpenFS = uicontrol( 'Parent', obj.controls.subframe2, 'Style', 'pushbutton', 'Position', [10 25 100 40], ...
                'String', 'Open', 'Callback', @obj.btnOpenFS, 'FontSize', 11 , 'FontWeight', 'bold');
            
            % text box inside subframe 2
            obj.controls.txtFS= uicontrol( 'Parent', obj.controls.subframe2, 'Style', 'edit', 'Position', [115 26 168 36], ...
                'FontSize', 12, 'string', {' FS ribbon (*.nii)'} , 'HorizontalAlignment', 'left', 'BackgroundColor', 'w' ,'enable','inactive');
            
            
            %% Inside subframe 3 - frame 1 =================================================
            % Button 1 inside subframe 3
            obj.controls.btnOpenCT1 = uicontrol( 'Parent', obj.controls.subframe3, 'Style', 'pushbutton', 'Position', [10 25 100 40], ...
                'String', 'Open', 'Callback', @obj.btnOpenCT1, 'FontSize', 11 , 'FontWeight', 'bold');
            
            % text box inside subframe 3
            obj.controls.txtCT1= uicontrol( 'Parent', obj.controls.subframe3, 'Style', 'edit', 'Position', [115 26 168 36], ...
                'FontSize', 12, 'string', {' CT scan (*.nii)'} , 'HorizontalAlignment', 'left', 'BackgroundColor', 'w' ,'enable','inactive');
            
            
            %% Inside subframe 4 - frame 2 ========================================================================================================================
            % Button 1 inside subframe 4
%             obj.controls.btnOpenCT2 = uicontrol( 'Parent', obj.controls.subframe4, 'Style', 'pushbutton', 'Position', [10 25 100 40], ...
%                 'String', 'Open', 'Callback', @obj.btnOpenCT1, 'FontSize', 11 , 'FontWeight', 'bold');
            
            % text box inside subframe 4
            obj.controls.txtCT2 = uicontrol( 'Parent', obj.controls.subframe4, 'Style', 'edit', 'Position', [10 26 168+105 36], ...
                'FontSize', 12, 'string', {' CT scan (*.nii)'} , 'HorizontalAlignment', 'left', 'BackgroundColor', 'w' ,'enable','inactive');
            
            %% Inside subframe 5 - frame 3
            %text box
            obj.controls.txtCV = uicontrol( 'Parent', obj.controls.subframe5, 'Style', 'text', 'Position', [10 145 213 25], ...
                'FontSize', 12, 'string', {'Electrode max. intensity (-5)'} , 'HorizontalAlignment', 'left','enable','inactive');
            
            %edit text box
            obj.controls.edtCV = uicontrol( 'Parent', obj.controls.subframe5, 'Style', 'edit', 'Position', [223 145 60 36], ...
                'FontSize', 12, 'string', {'3999'} ,'Callback', @obj.edtCV,  'HorizontalAlignment', 'center', 'BackgroundColor', 'w' ,'enable','on');
            
            %text box
            obj.controls.txtR = uicontrol( 'Parent', obj.controls.subframe5, 'Style', 'text', 'Position', [10 85 213 25], ...
                'FontSize', 12, 'string', {'Electrode volume (e.g., 3)'} , 'HorizontalAlignment', 'left','enable','inactive');
            
            %edit text box
            obj.controls.edtR = uicontrol( 'Parent', obj.controls.subframe5, 'Style', 'edit', 'Position', [223 85 60 36], ...
                'FontSize', 12, 'string', {'3'} , 'Callback', @obj.edtR, 'HorizontalAlignment', 'center', 'BackgroundColor', 'w' ,'enable','on');
            
            %text box
            obj.controls.txtIS = uicontrol( 'Parent', obj.controls.subframe5, 'Style', 'text', 'Position', [10 25 213 25], ...
                'FontSize', 12, 'string', {'Interelectrode space (e.g., 1)'} , 'HorizontalAlignment', 'left','enable','inactive');
            
            %edit text box
            obj.controls.edtIS = uicontrol( 'Parent', obj.controls.subframe5, 'Style', 'edit', 'Position', [223 25 60 36], ...
                'FontSize', 12, 'string', {'1'} , 'Callback', @obj.edtIS, 'HorizontalAlignment', 'center', 'BackgroundColor', 'w' ,'enable','on');
            
            
            %% Inside subframe 6 - frame 3
            %select a method
            
            obj.controls.Methods = uibuttongroup('Parent', obj.controls.subframe6, 'Visible', 'on', 'SelectionChangedFcn', @obj.radiobtnSelectionMethod, 'Position', [0 0.66 1 0.5], 'Bordertype', 'none');
            
            obj.controls.radiobtn0 = uicontrol( obj.controls.Methods, 'Style', 'radiobutton', 'Position', [20 75 230 25], ...
                'FontSize', 12, 'string', 'Trotta' ,'HandleVisibility','off','Value', 1, 'HorizontalAlignment', 'left','enable','on');

            
            
            
            %enter subject name:
            %text box
            obj.controls.txtSbjName = uicontrol( obj.controls.subframe6, 'Style', 'text', 'Position', [10 205 150 40], ...
                'FontSize', 12, 'string', {'Subject Name:'} , 'FontWeight', 'bold','HorizontalAlignment', 'left','enable','inactive');
            %edit text box
            obj.controls.edtSbjName =  uicontrol( 'Parent', obj.controls.subframe6, 'Style', 'edit', 'Position', [130 224 155 25], ...
                'FontSize', 12, 'string','NAME' , 'HorizontalAlignment', 'center', 'BackgroundColor', 'w' ,'enable','on');
            
            %select an hemisphere
            obj.controls.Hemisphere = uibuttongroup('Parent', obj.controls.subframe6, 'Visible', 'on', 'SelectionChangedFcn', @obj.radiobtnSelectionHemisphere,'Position', [0 0.33 1 0.33],'Bordertype', 'none');
            
            obj.controls.radiobtn3 = uicontrol( obj.controls.Hemisphere, 'Style', 'radiobutton', 'Position', [160 75 230 25], ...
                'FontSize', 12, 'string', 'Left' ,'HandleVisibility','off','Value', 1, 'HorizontalAlignment', 'left','enable','on');
            
            obj.controls.radiobtn4 = uicontrol( obj.controls.Hemisphere, 'Style', 'radiobutton', 'Position', [220 75 230 25], ...
                'FontSize', 12, 'string', 'Right' ,'HandleVisibility','off','Value', 0, 'HorizontalAlignment', 'left','enable','on');
            
            obj.controls.txtHemisphere = uicontrol( obj.controls.Hemisphere, 'Style', 'text', 'Position', [10 60 150 40], ...
                'FontSize', 12, 'string', {'Implanted hemisphere:'} , 'FontWeight', 'bold','HorizontalAlignment', 'left','enable','inactive');
            
  
            %% Logging frame
            obj.controls.logframe = uipanel( 'Parent', obj.mainFig, 'Units', 'pixels', 'Position', [10 10 980 frameheight-20+40], ...
                'Title', 'Action Log:', 'FontSize', 14, 'FontWeight', 'bold', 'BorderType', 'line', 'HighlightColor', [.5 .5 .5] );
            
            % Inside the Logging frame:
            % text box
            obj.controls.txtLog = uicontrol( 'Parent', obj.controls.logframe, 'Style', 'text','max',2, 'Position', [10 10 960 145+40], ...
                'FontSize', 12, 'string', {'> Welcome to ALICE!', '> Please create a directory to start ALICE from scratch, or locate the existing directory.'}, 'HorizontalAlignment', 'left', 'BackgroundColor', 'w' ,'enable','inactive');
            
        end
        
        function CreateDirectory( obj )
            
            obj.settings.originaldir = [pwd '/'];
            obj.settings.scriptspath  = [fileparts( mfilename('fullpath') ) '/'];
            
            if exist([pwd '/ALICE'])==7
                errordlg('A folder named ALICE already exists. Please choose another directory or rename the existing folder.');
                
            else
                %make panels visible
                set(obj.controls.frame1,'Visible', 'on');
                set(obj.controls.frame2,'Visible', 'on');
                set(obj.controls.frame3,'Visible', 'on');
                
                %create directories
                mkdir('ALICE');
                cd ./ALICE;
                obj.settings.currdir = [pwd '/'];
                mkdir('log_info');
                %data folder
                %mkdir('data'); Mike getting rid of this folder
                %cd ./data;
                mkdir('3Dclustering');
                mkdir('CT');
                mkdir('coregistration')
                mkdir('MRI');
                mkdir('FreeSurfer');
                cd(obj.settings.currdir);
                %locate matlab scripts
                addpath(genpath([ obj.settings.scriptspath 'MATLAB_scripts']));
                %locate afni scripts
                cd(obj.settings.currdir);
                obj.settings.scriptspath  = [fileparts( mfilename('fullpath') ) '/']; % (gets this .m file's parent directory)
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/alignCTtoT1_shft_res.csh'], [obj.settings.currdir 'coregistration/']);
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/3dclustering.csh'], [obj.settings.currdir '3Dclustering/']);
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/select_electrode.csh'], [obj.settings.currdir '3Dclustering/']);
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/open_afni_suma.csh'], [obj.settings.currdir '3Dclustering/']);
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/indexify_electrodes.csh'], [obj.settings.currdir '3Dclustering/']);
                
                cd(obj.settings.currdir);
                addpath(genpath(obj.settings.currdir));
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},sprintf('> ALICE successfully created (%s). Please proceed to Step 1.',obj.settings.currdir)});
                loggingActions(obj.settings.currdir,1,' > ALICE successfully created. Please proceed to Step 1.');
            end
        end
        
        function LocateDirectory( obj, folderName )
            
            if nargin < 2
                folderName = uigetdir('.', 'Please locate ALICE folder.');
            end
            
            if folderName~=0
                
                %make panels visible
                set(obj.controls.frame1,'Visible', 'on');
                set(obj.controls.frame2,'Visible', 'on');
                set(obj.controls.frame3,'Visible', 'on');
                
                obj.settings.currdir = [folderName '/'];
                cd(folderName);
                
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str(end-5:end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> ALICE directory located: ' obj.settings.currdir],'> Please proceed to Step 1.'});
                loggingActions(obj.settings.currdir,1,[' > ALICE directory located: ' obj.settings.currdir]);
                loggingActions(obj.settings.currdir,1,' > Please proceed to Step 1.');
                %if files in folder then read them:
                %MRI
                E = lsCell([obj.settings.currdir 'MRI/']);
                if length(E) > 2
                    N = E(3).name;
                    obj.settings.MRI = [obj.settings.currdir 'MRI/' N];
                    obj.settings.loaded(1) = 1;
                    set(obj.controls.txtMRI, 'string',['...' obj.settings.MRI(end-28:end)]);
                    %log
                    str = get(obj.controls.txtLog, 'string');
                    if length(str)>=obj.settings.NUM_LINES
                        str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                    end
                    set(obj.controls.txtLog, 'string',{str{:},['> MRI scan selected: ' obj.settings.MRI]});
                    loggingActions(obj.settings.currdir,1,[' > MRI scan selected: ' obj.settings.MRI]);
                end
                %FS
                E = [obj.settings.currdir 'FreeSurfer/t1_class.nii'];
                if exist(E)~=0
                    obj.settings.FS = [obj.settings.currdir 'FreeSurfer/t1_class.nii'];
                    obj.settings.loaded(2) = 1;
                    set(obj.controls.txtFS, 'string',['...' obj.settings.FS(end-28:end)]);
                    %log
                    str = get(obj.controls.txtLog, 'string');
                    if length(str)>=obj.settings.NUM_LINES
                        str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                    end
                    set(obj.controls.txtLog, 'string',{str{:},['> FS segmentation selected: ' obj.settings.FS]});
                    loggingActions(obj.settings.currdir,1,[' > FS segmentation selected: ' obj.settings.FS]);
                end
                %CT
                if isfield(obj.settings, 'CT')
                    E = obj.settings.CT;
                    if exist(E)~=0
                        obj.settings.loaded(3) = 1;
                        set(obj.controls.txtCT1, 'string',['...' obj.settings.CT(end-28:end)]);
                        set(obj.controls.txtCT2, 'string',['...' obj.settings.CT(end-46:end)]);
                        %log
                        str = get(obj.controls.txtLog, 'string');
                        if length(str)>=obj.settings.NUM_LINES
                            str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                        end
                        set(obj.controls.txtLog, 'string',{str{:},['> CT scan selected: ' obj.settings.CT]});
                        loggingActions(obj.settings.currdir,1,[' > CT scan selected: ' obj.settings.CT]);
                    end
                end
                
                if exist([obj.settings.currdir '/log_info/settings.mat'])==2
                    
                    load([obj.settings.currdir '/log_info/settings.mat']);
                    try
                        obj.settings.R = settings.R;
                        obj.settings.CV = settings.CV;
                        obj.settings.IS = settings.IS;
                        set(obj.controls.edtCV, 'String',obj.settings.CV);
                        set(obj.controls.edtR, 'String',obj.settings.R);
                        set(obj.controls.edtIS, 'String',obj.settings.IS);
                    end
                    
                    try
                        obj.settings.subject = settings.subject;
                        obj.settings.Hemisphere = settings.Hemisphere;
                        obj.settings.Method = settings.Method;
                        set(obj.controls.edtSbjName, 'String',obj.settings.subject);


                        if strcmp(obj.settings.Hemisphere, 'Left')
                            set(obj.controls.radiobtn3,'Value', 1);
                            set(obj.controls.radiobtn4,'Value', 0);
                        else
                            set(obj.controls.radiobtn3,'Value', 0);
                            set(obj.controls.radiobtn4,'Value', 1);
                        end
                    end
                    
                    %log
                    str = get(obj.controls.txtLog, 'string');
                    if length(str)>=obj.settings.NUM_LINES
                        str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                    end
                    set(obj.controls.txtLog, 'string',{str{:},['> Saved subject name, cluster settings and grid settings loaded.']});
                    loggingActions(obj.settings.currdir,1,[' > Saved subject name, cluster settings and grid settings loaded.']);
                end
                
                addpath(genpath(obj.settings.currdir));
                
            else %if no directory is selected
                disp('! WARNING: No directory selected.');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},'> WARNING: No directory selected.'});

            end
        end
        
        function AlignCTtoMRI( obj )
            
            cd(obj.settings.currdir);
            
            if obj.settings.loaded(1) == 1 && obj.settings.loaded(3) == 1
                %log before
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:}, '> Aligning CT to MRI. This procedure might take several minutes. Please wait...'});
                loggingActions(obj.settings.currdir,1,' > Aligning CT to MRI... This might take several minutes. Please wait...');
                pause(1);
                
                cd([obj.settings.currdir '/coregistration/']);
                
                %clean directory if already run before:
                if 0
                    system('rm *.nii');
                    system('rm *.1D');
                    system('rm temp_*');
                end
                
                %copy CT local
                copyfileSafe(obj.settings.CT, pwd);
                
                %select T1
                T1_path = obj.settings.MRI;

                %align:
                system(['tcsh -x alignCTtoT1_shft_res.csh -CT_path ' obj.settings.CT ' -T1_path ' T1_path]);
                loggingActions(obj.settings.currdir,1, [' > tcsh -x alignCTtoT1_shft_res.csh -CT_path ' obj.settings.CT ' -T1_path' T1_path]);
                cd(obj.settings.currdir);
                
                fd = fopen('./coregistration/status.txt');
                if fd > 0
                    S = fscanf(fd,'%s');
                else
                    warning('No coregistration/status.txt file found!');
                end
                
                if fd > 0 && strcmp(S(end-20:end),'alignmentsuccessfully')
                    
                    %log after
                    str = get(obj.controls.txtLog, 'string');
                    if length(str)>=obj.settings.NUM_LINES
                        str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                    end
                    set(obj.controls.txtLog, 'string',{str{:}, '> CT and MRI have been co-registered. Please check the results before proceeding.'});
                    loggingActions(obj.settings.currdir,1,' > CT and MRI have been co-registered. Please check the results before proceeding.');
                    
                    %display instruction box:
                    h = msgbox({['To check whether the alignment was successful,', ' please check the overlayed files in AFNI.'],...
                        ['For that use the buttons OVERLAY and UNDERLAY ', 'to specify the CT and MRI, respectively. ', ...
                        'Use the intensity button on the volume windows (9),', ' to decrease/increase the intensity of the overlay.'],...
                        [' '],['- If the alignment is good, please close AFNI and proceed to Step 2.'],...
                        ['- If the alignment is NOT good, please check the if the CT and MRI where correctly',...
                        ' converted to *.nii. Check orientation and voxel sizes too.']},'How to check if alignment was successful?', 'help');
                else
                    %display instruction box:
                    h = msgbox({['Alignment failed,', ' please check the input files.'],...
                        ['Please check the if the CT and MRI where correctly',...
                        ' converted to *.nii. Check orientation and voxel sizes too.']},'Alignment failed', 'error');
                end
                if fd > 0
                    fclose(fd);
                end
                
            else
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:}, '>! WARNING: CT and MRI scans are not loaded yet. Use the buttons in Step 1 to load them.'});
                loggingActions(obj.settings.currdir,1,' >! WARNING: CT and MRI scans are not loaded yet. Use the buttons in Step 1 to load them.');
            end
            
        end
        
        function ExtractClusters ( obj )
            
            cd(obj.settings.currdir);
            
            if obj.settings.loaded(3) == 1
                
                try
                    system(['rm ' obj.settings.currdir '/coregistration/temp_ANAT.nii']);
                catch
                end
                %system(['rm ' obj.settings.currdir '/coregistration/*CT*.nii']);
                
                dclust = fullfile(obj.settings.currdir, '3dclustering');
                
                %log before 3dclustering
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:}, '> Extracting electrode clusters. This procedure might take some minutes. Please wait...'});
                loggingActions(obj.settings.currdir,2,' > Extracting electrode clusters. This procedure might take some minutes. Please wait...');
                pause(1);
                
                %3dclustering
                cd([obj.settings.currdir '/3dclustering/']);
                
                if ~exist(fullfile(dclust, ['3dclusters_r' num2str(obj.settings.R) '_is' num2str(obj.settings.IS) '_thr' num2str(obj.settings.CV) '.nii']), 'file') && exist( obj.settings.CT ,'file')
                    
                    system(['tcsh -x 3dclustering.csh -CT_path ' obj.settings.CT ' -radius ' num2str(obj.settings.R) ' -interelectrode_space ' num2str(obj.settings.IS) ' -clip_value ' num2str(obj.settings.CV)]);
                    %log after 3dclustering
                    str = get(obj.controls.txtLog, 'string');
                    if length(str)>=obj.settings.NUM_LINES
                        str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                    end
                    set(obj.controls.txtLog, 'string',{str{:}, '> Electrode clusters extracted. Please check results and then close SUMA.'});
                    loggingActions(obj.settings.currdir,2,' > Electrode cluster extracted. Please check results and then close SUMA.');
                    system(['suma -i 3dclusters_r' num2str(obj.settings.R) '_is' num2str(obj.settings.IS) '_thr' num2str(obj.settings.CV) '.gii']);
                else
                    str = get(obj.controls.txtLog, 'string');
                    if length(str)>=obj.settings.NUM_LINES
                        str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                    end
                    set(obj.controls.txtLog, 'string',{str{:}, '>! ERROR: delete any 3dclusters_rX_isX_thrX.nii files in the /3Dclustering directory (this function cannot overwrite files) AND Check that ct_implant.nii is inside /CT folder.'});
                    loggingActions(obj.settings.currdir,2,' > ! ERROR: delete any 3dclusters_rX_isX_thrX.nii files in the /3Dclustering directory (this function cannot overwrite files) AND Check if ct_implant.nii is inside /CT folder.');
                end
                
                cd(obj.settings.currdir);
                
            else
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:}, '>! WARNING: CT was not loaded yet. Use the button in Step 2 to load it.'});
                loggingActions(obj.settings.currdir,2,' >! WARNING: CT was not loaded yet. Use the button in Step 2 to load it.');
                
                cd(obj.settings.currdir);
                
            end
            
        end
        
        function SelectElectrodes( obj )
            
            cd(obj.settings.currdir);
            
            if isfield(obj.settings, 'chans')
                chans = obj.settings.chans;
            else
                chans = [];
            end
            
            %check if file exists
            if exist([obj.settings.currdir '/3Dclustering/3dclusters_r' num2str(obj.settings.R) '_is' num2str(obj.settings.IS) '_thr' num2str(obj.settings.CV) '.nii']) ~=0
                
                %log before select electrodes
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:}, '> Time to select the electrodes! Please select them according to the leads order.'});
                loggingActions(obj.settings.currdir,2,' > Time to select the electrodes! Please select them according to the leads order.');
                
                %move to 3dclustering folder
                cd([obj.settings.currdir '/3Dclustering/']);
                
                %define input file
                clust_set  = ['3dclusters_r' num2str(obj.settings.R) '_is'...
                    num2str(obj.settings.IS) '_thr' num2str(obj.settings.CV) '.nii'];
                clust_surf = ['3dclusters_r' num2str(obj.settings.R) '_is' ...
                    num2str(obj.settings.IS) '_thr' num2str(obj.settings.CV) '.gii'];
                
                %open afni and suma:
                system(['tcsh -x open_afni_suma.csh -CT_path ' obj.settings.CT '  -clust_set ' clust_set ' -clust_surf ' clust_surf ], '-echo');
                
                
                obj.subGUI = selectElGUI( obj.settings, obj.controls, chans );
                
                %display instruction box:
                if 0
                    h = msgbox({['Time to select the electrodes!'], ['Please have the leads order (electrode layout) at hand.'],...
                        [' '],['AFNI and SUMA interfaces will open, ', 'and you will see the electrodes clusters both in the volume and the surface spaces.'], ...
                        [' '],['With right click you can select the electrode on the surface,', ' and with the left click on the volume space.'],...
                        ['Use the click+drag to navigate in SUMA, scroll to zoom in and out and scroll-lock for panning.'],...
                        [' '],['Please use the matlab pop-up window to select the electrodes.'],...
                        },'Time to select electrodes!', 'help');
                end
                
            else
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:}, ['>! WARNING: There is no 3dclusters_r' num2str(obj.settings.R) '_is' num2str(obj.settings.IS) '_thr' num2str(obj.settings.CV) '.nii file in the /3Dclustering folder. Please check the settings above or extract clusters first.']});
                loggingActions(obj.settings.currdir,2,[' >! WARNING: There is no 3dclusters_r' num2str(obj.settings.R) '_is' num2str(obj.settings.IS) '_thr' num2str(obj.settings.CV) '.nii file in the /3Dclustering folder. Please check the settings above or extract clusters first.']);
            end
        end
        
    end
    
    methods % for callbacks
        
        % Close request function for the main dialog.
        function figCloseRequest( obj, hObject, ~ )
            delete( hObject );
            %clc;
            disp('     ___    ');
            disp('    |   |   _ ');
            disp(' ___| __|__|_|__        ');
            disp('|    | |        |_          ');
            disp('|____| |________|_| |  Together everyone achieves more  |           ');
            disp(' ');
            
            
        end
        
        %locate directory
        function btnLocateDirectory( obj, hObject, ~ )
            
            obj.LocateDirectory
        end
        
        %create directory
        function btnCreateDirectory( obj, hObject, ~ )
            
            obj.CreateDirectory;
            
        end
        
        %open MRI
        function btnOpenMRI( obj, hObject, ~, pathname )
            if nargin < 4; [filename, pathname] = uigetfile('../*.nii'); else, filename=[]; end
            
            cd(obj.settings.currdir);
            
            %clean up previous MRI
            if 0 && length(dir([obj.settings.currdir '/MRI/']))>2
                choice = questdlg({[' '],['Other MRI scan(s) have have been found in ./MRI folder.'], ...
                    ['If you choose to delete file(s), the existing file(s) will be deleted and replaced by the new file!'],[' ']},...
                    'WARNING!', 'Delete old file(s)', 'Keep old file(s)', 'Delete old file(s)'); 
                
                switch choice
                    case 'Delete old file'
                        delete([obj.settings.currdir '/MRI/*.nii']);
                        %log
                        str = get(obj.controls.txtLog, 'string');
                        if length(str)>=obj.settings.NUM_LINES
                            str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                        end
                        set(obj.controls.txtLog, 'string',{str{:}, '>! WARNING: Deleting previously loaded MRI file.'});
                end
            end
            
            
            
            if pathname~=0
                copyfileSafe([pathname filename], [obj.settings.currdir 'MRI/' filename ]);
                if isempty(filename)
                    [~, filename, ext] = fileparts(pathname);
                    filename = [filename ext];
                end
                obj.settings.MRI = [obj.settings.currdir 'MRI/' filename];
                obj.settings.loaded(1) = 1;
                settings = obj.settings; %#ok<NASGU>
                save([obj.settings.currdir '/log_info/settings'], 'settings');
                
                set(obj.controls.txtMRI, 'string',['...' obj.settings.MRI(end-28:end)]);
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> MRI scan selected: ' obj.settings.MRI]});
                loggingActions(obj.settings.currdir,1,[' > MRI scan selected: ' obj.settings.MRI]);
                
            else
                disp('! WARNING: MRI scan not selected.');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:}, '>! WARNING: MRI scan not selected.'});
                loggingActions(obj.settings.currdir,1,' >! WARNING: MRI scan not selected.');
            end
        end
        
        %open FS
        function btnOpenFS( obj, hObject, ~ , PathName)
            if nargin < 4; [FileName, PathName] = uigetfile('../*.nii'); else, FileName=[]; end
            
            cd(obj.settings.currdir);
            

            if PathName~=0
                %rename FS ribbon to t1_class.nii
                copyfileSafe([PathName FileName], [obj.settings.currdir 'FreeSurfer/t1_class.nii']);
                obj.settings.FS = [obj.settings.currdir 'FreeSurfer/t1_class.nii'];
                obj.settings.loaded(2) = 1;
                settings = obj.settings;
                save([obj.settings.currdir '/log_info/settings'], 'settings');
                
                set(obj.controls.txtFS, 'string', ['...' obj.settings.FS(end-28:end)]);
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> FS segmentation selected: ' obj.settings.FS]});
                loggingActions(obj.settings.currdir,1,[' > FS segmentation selected: ' obj.settings.FS]);
            else
                disp('! WARNING: FreeSurfer segmentation file not selected.');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},'>! WARNING: FreeSurfer segmentation file not selected.'});
                loggingActions(obj.settings.currdir,1,' >! WARNING: FreeSurfer segmentation file not selected.');
            end
        end
        
        %Open CT
        function btnOpenCT1( obj, hObject, ~, PathName )
            if nargin < 4; [FileName, PathName] = uigetfile('../*.nii'); else, FileName=[]; end
            
            
            cd(obj.settings.currdir);
            
            
            
            
            if PathName~=0
                %rename CT to CT_highresRAI.nii
                if ~strcmpi([PathName FileName], fullfile(obj.settings.currdir, 'CT', FileName))
                    copyfileSafe([PathName FileName], fullfile(obj.settings.currdir, 'CT'));
                end
                
                if isempty(FileName)
                    [~,FileName] = fileparts(PathName);
                    FileName = [FileName '.nii'];
                end
            
                fname_ct = fullfile(obj.settings.currdir, 'CT', FileName);
                obj.settings.CT = fname_ct;
                obj.settings.loaded(3) = 1;
                
                %update text boxes
                set(obj.controls.txtCT1, 'string',['...' obj.settings.CT(end-28:end)]);
                set(obj.controls.txtCT2, 'string', ['...' obj.settings.CT(end-34:end)]);
                
                %extract ct max value
                system(['3dBrickStat -slow ' fname_ct ' > temp_ct_val.txt']);
                Faux = fopen('temp_ct_val.txt');
                obj.settings.CV = fscanf(Faux, '%d')-5;
                set(obj.controls.edtCV, 'string', {num2str(obj.settings.CV)});
                fclose(Faux);
                delete(['temp_ct_val.txt']);
                
                %save settings
                settings = obj.settings;
                save([obj.settings.currdir '/log_info/settings'], 'settings');
                
                %log CT
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> CT scan selected: ' obj.settings.CT]});
                loggingActions(obj.settings.currdir,1,[' > CT scan selected: ' obj.settings.CT]);
                loggingActions(obj.settings.currdir,2,[' > CT scan selected: ' obj.settings.CT]);
                
                %log CV
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> Electrode maximum intensity selected: ' num2str(obj.settings.CV)]});
                loggingActions(obj.settings.currdir,2,[' > Electrode maximum intensity selected: ' num2str(obj.settings.CV)]);
                
            else
                disp('! WARNING: CT scan not selected.');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},'>! WARNING: CT scan not selected.'});
                loggingActions(obj.settings.currdir,1,' >! WARNING: CT scan not selected.');
                loggingActions(obj.settings.currdir,2,' >! WARNING: CT scan not selected.');
            end
        end
        
        %Align
        function btnAlignCTtoMRI( obj, hObject, ~ )
            
            ojb.controls.btnAlignCTtoMRI.Enable = 'off';
            obj.AlignCTtoMRI;
            ojb.controls.btnAlignCTtoMRI.Enable = 'on';
            
        end
        
        %enter clip value
        function edtCV( obj, hObject, ~ )
            
            str = get(obj.controls.edtCV, 'string');
            obj.settings.CV = str2num(str{1});
            
            if isempty(obj.settings.CV)
                disp('! WARNING: Maximum intensity must be an integer value bigger than 0.');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},'>! WARNING: Maximum intensity must be an integer value bigger than 0.'});
                loggingActions(obj.settings.currdir,2,' >! WARNING: Maximum intensity must be an integer value bigger than 0.');
            else
                settings = obj.settings;
                save([obj.settings.currdir '/log_info/settings'], 'settings');
                
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> Electrode maximum intensity selected: ' num2str(obj.settings.CV)]});
                loggingActions(obj.settings.currdir,2,[' > Electrode maximum intensity selected: ' num2str(obj.settings.CV)]);
            end
            
        end
        
        %enter radius
        function edtR( obj, hObject, ~ )
            
            str = get(obj.controls.edtR, 'string');
            obj.settings.R = str2double(char(str));
            
            if isempty(obj.settings.R)
                disp('! WARNING: Electrode volume must be an integer value bigger than 0.');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},'>! WARNING: Electrode volume must be an integer value bigger than 0.'});
                loggingActions(obj.settings.currdir,2,' >! WARNING: Electrode volume must be an integer value bigger than 0.');
            else
                settings = obj.settings;
                save([obj.settings.currdir '/log_info/settings'], 'settings');
                
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> Electrode volume selected: ' num2str(obj.settings.R)]});
                loggingActions(obj.settings.currdir,2,[' > Electrode volume selected: ' num2str(obj.settings.R)]);
            end
        end
        
        %enter interelectrode distance
        function edtIS( obj, hObject, ~ )
            
            str = get(obj.controls.edtIS, 'string');
            obj.settings.IS = str2double(char(str));
            
            if isempty(obj.settings.IS)
                disp('! WARNING: Interelectrode space must be an integer value bigger than 0.');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},'>! WARNING: Interelectrode space must be an integer value bigger than 0.'});
                loggingActions(obj.settings.currdir,2,' >! WARNING: Interelectrode space must be an integer value bigger than 0.');
            else
                settings = obj.settings;
                save([obj.settings.currdir '/log_info/settings'], 'settings');
                
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> Interelectrode space selected: ' num2str(obj.settings.IS)]});
                loggingActions(obj.settings.currdir,2,[' > Interelectrode space selected: ' num2str(obj.settings.IS)]);
            end
            
        end
        
        %Select electrodes
        function btnExtractClusters( obj, hObject, ~ )
            
            
            ojb.controls.btnExtractClusters.Enable = 'off';
            obj.ExtractClusters;
            ojb.controls.btnExtractClusters.Enable = 'on';
            
        end
        
        %Select electrodes
        function btnSelectElectrodes( obj, hObject, ~ )
            
            obj.controls.btnSelectElectrodes.Enable = 'off';
            obj.SelectElectrodes;
            obj.controls.btnSelectElectrodes.Enable = 'on';
            
        end
        
        
        %choose hemisphere
        function radiobtnSelectionHemisphere( obj, hObject, callbackdata )
            
            obj.settings.Hemisphere = callbackdata.NewValue.String;
            
            %log
            str = get(obj.controls.txtLog, 'string');
            if length(str)>=obj.settings.NUM_LINES
                str = str( (end - (obj.settings.NUM_LINES-1)) :end);
            end
            set(obj.controls.txtLog, 'string',{str{:},['> ' obj.settings.Hemisphere ' hemisphere selected.']});
            loggingActions(obj.settings.currdir,3,[' > ' obj.settings.Hemisphere ' hemisphere selected.']);
            
        end
        
      
        
        %visualize results
        function btnVisualize( obj, hObject, ~ )
            
            cd(obj.settings.currdir);
            
            %if no FS
            if ~isfield(obj.settings, 'FS')
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},'>! WARNING: Please load a FreeSurfer segmentation!'});
                loggingActions(obj.settings.currdir,3,' >! WARNING: Please load a FreeSurfer segmentation!');
                return;
            end
            
            subject = get(obj.controls.edtSbjName, 'String');
            
            %if subject name is empty
            if isempty(subject) || strcmp(subject,'')
                obj.settings.subject = ' ';
                disp('>! WARNING: No name entered.');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},'>! WARNING: No name entered. Try again!'});
                loggingActions(obj.settings.currdir,3,' >! WARNING: No name entered. Try again!');
                status = 0;
                return;
            
            else
                obj.settings.subject = subject;
                settings = obj.settings;
                save([obj.settings.currdir '/log_info/settings'], 'settings');
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> Name entered: ' subject]});
                loggingActions(obj.settings.currdir,3, [ ' > Name entered: ' subject]);
            end
            
            if strcmp(obj.settings.Method, 'Trotta')
                disp(['> Applying ' obj.settings.Method '... Please wait until a figure with the projected electrodes appears.']);
                %log
                str = get(obj.controls.txtLog, 'string');
                if length(str)>=obj.settings.NUM_LINES
                    str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                end
                set(obj.controls.txtLog, 'string',{str{:},['> Applying ' obj.settings.Method '... Please wait until a figure with the projected electrodes appears.']});
                loggingActions(obj.settings.currdir,3,[' > Applying ' obj.settings.Method '... Please wait until a figure with the projected electrodes appears.']);
                pause(1);
                
                % run_trotta_projection
                status = 0;
                
                if status==1
                    disp('> Electrode projection completed. Please find the results in ./results_trotta/projected_electrodes_coord/.');
                    %log
                    str = get(obj.controls.txtLog, 'string');
                    if length(str)>=obj.settings.NUM_LINES
                        str = str( (end - (obj.settings.NUM_LINES-1)) :end);
                    end
                    set(obj.controls.txtLog, 'string',{str{:},'> Electrode projection completed. Please find the results in ./results_trotta/projected_electrodes_coord/.'});
                    loggingActions(obj.settings.currdir,3,' > Electrode projection completed. Please find the results in ./results_trotta/projected_electrodes_coord/.');
                    
                end
            end
            
        end
        
        
    end
end

function varargout = copyfileSafe(varargin)
    % copyfileSafeSafe uses a unix command to copy the file iff a preliminary matlab copyfileSafe fails
    %
    % Usage: [SUCCESS, MESSAGE] = copyfileSafeSAFE(SOURCE,DESTINATION,MODE)
    %
    % See Also: copyfileSafe
    
    [success, message] = copyfile(varargin{:});
    if ~success
        [status, result] = unix(sprintf('cp %s %s', varargin{1}, varargin{2}));
        success = (status == 0);
        message = result;        
    end
        
    varargout = {success, message};
end