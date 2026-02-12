% script localizer_rois_fsAve
%
%   run this script (re)generate the ROI_mesh_LUT for average brains (fsAve or colins_N27)
%   
%
%   File Output:
%       eeg_toolbox/visualize/data/fsaverage|colins_N27/ROIC_mesh_LUT_lh/rh.csv
%
%
%   6/2020 created by JW (by chopping down localizer_rois) to increase the ROIC to mesh radius for fsAve 
%



    %TOOLBOX        = '/Users/youssefdm/Documents/eeg_toolbox/';        % eeg_toolbox directory
    TOOLBOX         = getToolboxDir;                                    %  JW changed to function call so anybody can use this
    N_ROIC          = 2400;                                             % (coarse=600, fine=2400)
    R_ROI_MAX_DIST  = 30;                                               % maximum distance to store to vertices from ROIC's (was 10, JW changed to 30)

    
    
    WHICH_BRAIN = 2; %- 1 for fsAve,  2 for N27
    
    
    if WHICH_BRAIN==1,
        % fsaverage
        subj = 'fsaverage';
        locDirs.fs_subj = fullfile(TOOLBOX,'visualize/data/',subj);
        locDirs.roi     = locDirs.fs_subj;
        sumaDir         = 'SUMA';
        
    else
        % colin
        subj = 'colin_N27';
        locDirs.fs_subj = fullfile(TOOLBOX,'visualize/data/',subj);
        locDirs.roi     = locDirs.fs_subj;
        sumaDir         = 'suma_MNI';
        
    end
    
    
    
    hems = {'rh' 'lh'}; 
    hems=hems(~cellfun('isempty',hems));
    hems=char(hems);
    hems=hems(:,1)'; % char array of r/l or rl
    
    
    
    for hemisphere = hems
        fname_roi_mesh_lut_csv = fullfile(locDirs.roi, sprintf('ROIC_mesh_LUT_%sh.csv', hemisphere));
        fname_roi_mesh_lut_mat = fullfile(locDirs.roi, sprintf('ROIC_mesh_LUT_%sh.mat', hemisphere));
        % This table contains a VERY LARGE lookup table used to map ROICs to their region with columns:
        % 
        % note for Mike: the table is ROIC --> mesh. We actually only need electrode --> ROIC --> mesh
        %                ie it can just use union(ROIC_verts) as a starting set
        %
        % CSV table columns:
        %     vertex       :     mesh index within MAX_RADIUS of an ROIC mesh
        %     d            :     distance from mesh intext to ROIC mesh index
        %     ROIC_mesh_ndx:     ROIC mesh index under comparison
        %     ROIC_ndx     :     ROIC index number (from 1 to 2400) under evaluation
        % "ROIC_mesh_LUT_lh.csv (and then rh version): contents are ROIC index / ROIC vertex index / neighboring vertex index / neighbor distance"
        %
        % USAGE: braindata2 loads into roic_roi_lh/rh, and this table is used to map verticies to ROICs (vertex2ROI; for aggregating effects across subj), and ROICs to verticies (roic2roi; for visualizing)
        % USAGE: ES and JW used this to create ROIC-to-ROIC adjacency matrix of fsAverage for cluster correction, and to quantify ROIC to ROIC distance in ROIeffectFromElectrodes.m (ES custom code)
        
        
        if ~exist(fname_roi_mesh_lut_mat, 'file') 

            fname_default_ROIC  = fullfile(TOOLBOX, sprintf('visualize/data/zLocalize_ROI_centers_%sh.mat', hemisphere));
            pial_filename       = fullfile(locDirs.fs_subj,sumaDir,['std.141.' hemisphere 'h.pial.gii']);
            assert(exist(pial_filename, 'file') > 0, 'File not found: %s\n', pial_filename);

            pial_surf           = gifti(pial_filename);
            

            % ROIC standard file (There is 1 global standard per hemisphere; ie subject-independent)
            ROIC_filename   = uigetfile_if_not_exist(fname_default_ROIC, '*.mat', 'Select ROI Centers file');
            try
                temp        = load(ROIC_filename);
                std_ROICs   = temp.vert_idx(1:N_ROIC);
                clear('temp');
            catch e
                fprintf('Error: could not load ROICs: %s\n', e.message);
            end

            %- do this for all ROICs (for standard brain)
            %relevant_ROICs = sort(std_ROICs); %- dont sort... screws things up!
            relevant_ROICs = std_ROICs;
            
            
            % ROIC --> ROI growth
            fprintf('Calculating ROIC->ROIs %sh...\n', hemisphere);
            % create ROIC_mesh_LUT_%sh.csv / braindata.roi.roic_roi_%sh / fname_roi_mesh_lut:
            roi_mesh_d_lut_hem = grow_ROIs(pial_surf, relevant_ROICs, R_ROI_MAX_DIST); 

            %writetable(roi_mesh_d_lut_hem, fname_roi_mesh_lut_csv); % per hemisphere... will eventually cut this. no reason do use csv
            save(fname_roi_mesh_lut_mat, 'roi_mesh_d_lut_hem','-v7'); % per hemisphere... about 1/4 the space on disk as csv
       
        else
            fprintf('\n %s already exists.  \nIf you want to make a new one, delete or move the old one and rerun.',fname_roi_mesh_lut);
        end 
        

    
    end

    