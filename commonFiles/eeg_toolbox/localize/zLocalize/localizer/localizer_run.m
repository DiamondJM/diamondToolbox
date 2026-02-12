% localizer_run
%
%
% Revision History
%   10/2017 MST - Created. Some code is built-in for the old_tal-->new_tal conversion
%   6/18/2020 SJ: changed readtables to readtableSafe to be compatible with MATLAB 2020a
%                 - implement fix for all_coords and all_dural_coords
%                   variables- for some reason horzcat is not working with
%                   2020a

% Dev Notes
%  - Look into dcm2niix_afni for a MRI_raw --> dicom/nifti option
%  - Look into hippocampal subfields freesurfer package
%  - switch to fat_proc_convert_dcm_anat? (see https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/tutorials/fatcat_prep/main_toc.html)
%  - use a fat_proc to axialize a T2 (if we have a T2 Coronal)
%  - look into this: "Using nifti_tool to remove ID info, if possible" for anonymization lead
%
function localizer_run(subj, rootEEGdir, numCTs)
    directory = evalin('base','dir');
    if nargin < 1, subj = 'NIH046'; end
    if nargin < 2, rootEEGdir = directory.loctest; end
    if nargin < 3, numCTs = 1; end

    
    % ---------- Settings ---------- %
    USE_SYM_LINK    = 1; 
    FS_BIN          = '/Applications/freesurfer/7.1.0/bin/';            %SJ      % '/Applications/freesurfer/bin/'
    AFNI_BIN        = '/Users/jacksonsn/abin';                           % '/Users/trottams/abin'; '/Users/youssefdm/abin/'
    SLICER_BINARY   = '/Applications/Slicer.app/Contents/MacOS/Slicer'; % '/Applications/Slicer.app/Contents/MacOS/Slicer'
    DURAL_SNAP      = 0;                                                % 1 to use snap-to-dura post-projection as final. 0 to record in separate snap file
    USE_T2          = 0;                                                % NOT YET IMPLEMENTED
    BATCH           = 0;                                                % If 1, don't prompt
    SHOW_FIGURES    = 1;                                                % whether to launch figures
    
    USE_ALICE       = 1;
    RERUN_SURF      = 1;
    RERUN_SUMA      = 0;
    RERUN_ANCHORS   = 0;
    RERUN_PROJ      = 0;
    RERUN_XFM       = 0;
    ONLY_CT_NUM     = 0;	% 0 default. 1 or 2 for just a particular CT
    
    addpathSystem('/Users/jacksonsn/abin'); %'/Users/youssefdm/abin'
    
    RERUN_ROI       = 0;
    RERUN_ATLAS     = 1;
    fprintf('\n+--------------------------------------------------------------------------------------------------+\n');
    fprintf('|                                  ** zLocalize Pipeline **                                        |\n');
    fprintf('|                                      localizer_run.m                                             |\n');
    fprintf('+--------------------------------------------------------------------------------------------------+\n');
    
    %fprintf('---------------------------------------------------------------------------------------------------\n');
    fprintf('|Parameters \n');
    fprintf('|\t %12s:\t %s\n', 'Subject', subj);
    fprintf('|\t %12s:\t %s\n', 'rootEEGdir', rootEEGdir);
    fprintf('|\t %12s:\t %d\n', 'numCTs', numCTs);
    fprintf('|\t %12s:\t %d\n', 'USE_ALICE', USE_ALICE);
    fprintf('|\t %12s:\t %d\n', 'RERUN_SURF', RERUN_SURF);
    fprintf('|\t %12s:\t %d\n', 'RERUN_SUMA', RERUN_SUMA);
    fprintf('|\t %12s:\t %d\n', 'RERUN_ANCHORS', RERUN_ANCHORS);
    fprintf('|\t %12s:\t %d\n', 'RERUN_PROJ', RERUN_PROJ);
    fprintf('|\t %12s:\t %d\n', 'RERUN_XFM', RERUN_XFM);
    fprintf('|\t %12s:\t %d\n', 'RERUN_ROI', RERUN_ROI);
    fprintf('|\t %12s:\t %d\n', 'RERUN_ATLAS', RERUN_ATLAS);
    fprintf('|\t %12s:\t %d\n', 'DURAL_SNAP', DURAL_SNAP);
    fprintf('|\t %12s:\t %d\n', 'BATCH', BATCH);
    fprintf('|\t %12s:\t %d\n', 'SHOW_FIGURES', SHOW_FIGURES);
    fprintf('|\t %12s:\t %d\n', 'ONLY_CT_NUM', ONLY_CT_NUM);
    
    %fprintf('+--------------------------------------------------------------------------------------------------\n');
    if ~BATCH
        if ~inputYN('Continue?')
            return;
        end
    end
    
    % ----------------------------------------------
    %% ---------- Create base directories ----------
    % ----------------------------------------------
    tal = fullfile(rootEEGdir, subj, 'tal');
    if ~exist(tal,'dir'), mkdir(tal); end
    zloc = fullfile(tal, 'zloc');
    if ~exist(zloc,'dir'), mkdir(zloc); end
    locDirs = localizer_create_directories(subj, rootEEGdir, numCTs);

    % ----------------------------------------------------
    %% ---------- Define variables (like paths) ----------
    % ----------------------------------------------------
    fname_default_ct_coord      = fullfile(locDirs.ct_1, 'coords.csv');
    fname_default_fs_subj       = fullfile(locDirs.fs, subj);
    fname_default_mr_nii        = fullfile(locDirs.mr_pre, 'mr_pre.nii');
    fname_default_mr_nii_t2     = fullfile(locDirs.mr_t2, 'mr_pre_t2.nii');
    fname_default_ct_nii        = fullfile(locDirs.ct_1, 'ct_implant.nii');
    fname_default_ct_nii_2      = fullfile(locDirs.ct_2, 'ct_implant.nii');
    element_info                = getElementInfo(subj, rootEEGdir);
    jacktable                   = getJackTable(subj, rootEEGdir);
    hems                        = [];

    % -------------------------------------------
    %% ---------- File lookup/creation ----------
    % -------------------------------------------
    
    % copy essential files from rootEEGdir_ORIG to rootEEGdir
    if 0 && COPY_ORIG
        copyOrigs(rootEEGdir, locDirs, subj, TAL_ORIG, rootEEGdir_ORIG, USE_ALICE, numCTs, TAL_ORIG_SUB)
    end 
        
    has_mr_nii  = exist(fname_default_mr_nii, 'file');
    has_mr_nii_t2 = exist(fname_default_mr_nii_t2, 'file');
    has_ct_nii = exist(fname_default_ct_nii, 'file');
    has_ct_nii_2 = exist(fname_default_ct_nii_2, 'file');
    mrT1_path = [];
    mrT2_path = [];
    ct1_path = [];
    ct2_path = [];
    
    % locate files
    fs_subj_dir = uigetfile_if_not_exist(fname_default_fs_subj, [], 'Choose freesurfer subject directory', USE_SYM_LINK, 1, BATCH);
    if ~has_mr_nii
        mrT1_path = uigetfile_if_not_exist(locDirs.mr_pre_dicom, [], 'Choose MR T1 dicom directory', 0, 1, BATCH);
    end
    if ~has_mr_nii_t2 && USE_T2
        mrT2_path = uigetfile_if_not_exist(locDirs.mr_t2_dicom, [], 'Choose MR T2 dicom directory', 0, 1, BATCH);
    end
    if ~has_ct_nii
        ct1_path = uigetfile_if_not_exist(locDirs.ct_1_dicom, [], 'Choose 1st CT dicom directory', 0, 1, BATCH);
    end
    if numCTs == 2 && ~has_ct_nii_2
        ct2_path = uigetfile_if_not_exist(locDirs.ct_2_dicom, [], 'Choose 2nd CT dicom directory', 0, 1, BATCH);
    end
    
    if isempty(mrT1_path) || is_empty_dir(mrT1_path)
        uigetfile_if_not_exist(fname_default_mr_nii, '*.nii', 'Select MR nifti', 0, 0, BATCH);
        has_mr_nii_t2 = exist(fname_default_mr_nii, 'file');
    end
    
    if USE_T2 && (isempty(mrT2_path) || is_empty_dir(mrT2_path))
        uigetfile_if_not_exist(fname_default_mr_nii_t2, '*.nii', 'Select MR T2 nifti', 0, 0, BATCH);
    end
    
    if isempty(ct1_path) || is_empty_dir(ct1_path)
        uigetfile_if_not_exist(fname_default_ct_nii, '*.nii', 'Select CT nifti', 0, 0, BATCH);
    end
    if numCTs == 2 && ( isempty(ct2_path) || isempty(ct1_path) )
        uigetfile_if_not_exist(fname_default_ct_nii, '*.nii', 'Select CT2 nifti', 0, 0, BATCH);
    end

    
    % display located files
    localizer_get_status(subj, rootEEGdir, 1);
    
    % -----------------------------------------------------
    %% ---------- Process the DICOMs into NIFTIs ----------
    % -----------------------------------------------------
    if ~has_mr_nii
        read_dicoms(locDirs.mr_pre_dicom, 'mr_pre',locDirs.mr_pre, 'afni_bin',AFNI_BIN, 'freesurfer_bin',FS_BIN);
    end
    if ~has_mr_nii_t2 && USE_T2
         read_dicoms(locDirs.mr_t2_dicom, 'mr_pre',locDirs.mr_t2, 'afni_bin',AFNI_BIN, 'freesurfer_bin',FS_BIN);
    end
    if ~has_ct_nii
        fprintf('Looking for CT dicoms in %s...', locDirs.ct_1_dicom);
        if is_empty_dir(locDirs.ct_1_dicom)
            fprintf('No Dicoms found.\n');
        else
            read_dicoms(locDirs.ct_1_dicom, 'ct_implant', locDirs.ct_1); % this will have to be done once for every CT
            fprintf('done!\n');
        end
    end
    if numCTs == 2 && ~exist(fullfile(locDirs.ct_2, 'ct_implant.nii'), 'file')
        read_dicoms(locDirs.ct_2_dicom, 'ct_implant', locDirs.ct_2); % this will have to be done once for every CT
    end
    


    % -----------------------------------------
    %% ---------- Create Surfaces -------------
    % -----------------------------------------
    status = localizer_get_status(subj, rootEEGdir, 0);
    fname_t2 = [];
    if USE_T2
        if status.n_mr_nii2 > 0
            fname_t2 = fullfile(locDirs.mr_t2, 'mr_pre_t2.nii');
            assert(exist(fname_t2, 'file') > 0, 'File not found (make sure it is named exactly this): %s', fname_t2);
        else
            warning('Cannot find T2, continuing without using T2');
        end
    end
        
    
    if numel(lsCell(locDirs.fs_subj)) < 1 || RERUN_SURF
        create_surf(subj,fullfile(locDirs.mr_pre,'mr_pre.nii'), locDirs.fs, 'freesurfer_bin',FS_BIN, 'T2_nifti', fname_t2);
        % Did this step run successfully and *completely* ?
        % => you should have these files in zloc/freesurfer/[subj]/surf/ :
        %   -  lh.pial-outer-smoothed
        %   -  rh.pial-outer-smoothed
        % They must be present before calling suma(...)
    else
        fprintf('Subject freesurfer directory non-empty; skipping surface reconstruction.\n');
    end

    % this will not do much if they already exist
    if numel(lsCell(fullfile(locDirs.fs_subj, 'SUMA'))) < 1 || RERUN_SUMA
        
        suma(subj, locDirs.fs, 'afni_bin',AFNI_BIN, 'freesurfer_bin',FS_BIN, 'rerun',RERUN_SUMA);
        % Did this step run successfully and *completely* ?
        % => you should have these files in zloc/freesurfer/[subj]/SUMA/ :
        %   - rh.pial-outer-smoothed.gii
        %   - lh.pial-outer-smoothed.gii
    else
        fprintf('Subject SUMA directory non-empty; skipping SUMA step.\n');
    end
        
    % read hemispheres from element_info
    if ~isempty(element_info)    
        temp = char(unique(element_info.whichHemi(~cellfun(@isempty, element_info.whichHemi))));
        if isempty(temp)
            fprintf('\nError: Could not determine the hemisphere from element_info. Aborting now.\n');
            return;
        else
            hems=temp(:,1)'; % first character: l/r
        end
    else
        hems = [];
    end
    
    for ctNum = 1:numCTs
 
        if ONLY_CT_NUM > 0 && ctNum ~= ONLY_CT_NUM
            continue;
        end
        
        fprintf('\n----------------------------------------------------------------------------------------------------\n');
        fprintf('Running CT %d of %d\n', ctNum, numCTs);
        fprintf('----------------------------------------------------------------------------------------------------\n');
        
        % check if we have projected xyz for every lead already (we assume if we find all for CT2, we also found them for CT1)
        run_projection  = 1;
        %chans2project   = getLeads(subj, rootEEGdir, 'hardwareType',{'subdural','micro-subdural','depth'}, 'implant',numCTs, 'chanType','PHYS' , 'jacktable',jacktable);
        fname           = fullfile(locDirs.mr_pre, sprintf('coords_%d.csv', ctNum));
        if exist(fname, 'file')
            fprintf('Skipping projection -- found file: %s\n', fname);
            run_projection = 0;
        end
        if ~run_projection && ~RERUN_PROJ && ~RERUN_ANCHORS
            continue;
        end
        
        % check if we have an unprojected xyz for every lead already
        run_xfm = 1;
        ct_leads = getLeads(subj, rootEEGdir, 'implant', ctNum, 'jacktable',jacktable);
        ct_dir = locDirs.(sprintf('ct_%d_xfm',ctNum));
        if ctNum == 1
            ct2_str = '';
        else
            ct2_str = '_2';
        end
        fname = fullfile(ct_dir, 'coords.csv');
        if exist(fname, 'file')
            leads = readtableSafe(fname); %SJ - changed from readtable
            if numel(intersect(leads.chanName, ct_leads)) == numel(ct_leads) || ...
                    ctNum == 2 && ~isempty(leads)
                fprintf('Skipping transform: Found all leads for this CT marked in %s\n', fname);
                run_xfm = 0;
            else
                fprintf('The following leads were not found in the transformed-to-MR coords (%s):\n', fname);
                disp(setdiff(ct_leads, leads.chanName));
                fprintf('This most-likely means that you forgot to mark these leads. It could mean these leads are sh channels or \n');
                fprintf('should be marked recDontUse. In special cases, you can allow the exception by breaking here to manually set run_xfm=0.\n');
                fprintf('Consider using the function copyInterhemXYZ.m for interhemispheric elextrodes.\n');
                fprintf('By default, this pipeline will proceed to launch ALICE so that you can select more leads.\n');
                fprintf('Break NOW if you don''t want to continue...');
                for t=10:-1:1, fprintf('%d ',t); pause(1); end
                disp('');
            end
        end
        
        
        
        
        if USE_ALICE
            % ----------------------------------------------------------------------
            %% ---------- Run ALICE (GUI for CT selection + MR transform) ----------
            % ----------------------------------------------------------------------
            if run_xfm || RERUN_XFM
                for hemisphere = hems
                    runAgain = 1;
                    while runAgain
                        runAgain = 0;
                        try
                            run_alice(zloc,subj, locDirs, ctNum, hemisphere)
                        catch e
                            fprintf('run_alice error: %s (line %d)\n', e.message, e.stack.line);
                            fprintf('(Note: run_alice is called in a try-catch loop. Break before its call and execute manually to investigate errors)\n');
                            runAgain = inputYN(sprintf('Try to run again for hemisphere %s? ', hemisphere));
                        end
                    end
                end

                % read the MR-coordinate text files and convert to CSV. Make sure to use the most recent entries
                if numel(hems) == 1
                    fname = fullfile(ct_dir, sprintf('coords_%s.txt', hems(1)));
                    fd = fopen(fname);
                    assert(fd > 0, 'Error opening %s', fname);
                    x = textscan(fd, '%f %f %f %s\n');
                    fclose(fd);
                    t = table(x{[4,1:3]}, 'variableNames', {'chanName', 'x','y','z'});
                    [~,ndx] = unique(flip(t.chanName), 'stable'); ndx = height(t)-ndx+1; % flip=most recent
                    t = t(ndx, :);
                    writetable(t, fullfile(ct_dir, 'coords.csv'));
                else
                    t = table;
                    chanName_all = [];
                    xyz_all = [];
                    for i = 1:2
                        fname = fullfile(ct_dir, sprintf('coords_%s.txt', hems(i)));
                        fd = fopen(fname, 'r');
                        if fd <= 0
                            warning('Could not find %s', fname);
                            continue;
                        end
                        x = textscan(fd, '%f %f %f %s\n');
                        fclose(fd);
                        chanName = x{4};
                        xyz = cat(2,x{1:3});
                        [~,ndx] = unique(flip(chanName), 'stable'); ndx=numel(chanName)-ndx+1;% flip = most recent
                        chanName = chanName(ndx);
                        xyz = xyz(ndx,:);
                        chanName_all = cat(1, chanName_all, chanName);
                        xyz_all = cat(1, xyz_all, xyz);
                    end
                    if isempty(xyz_all)
                        warning('No electrodes marked');
                    else
                        t = table(chanName_all, xyz_all(:,1), xyz_all(:,2), xyz_all(:,3), 'variableNames', {'chanName', 'x','y','z'});
                        filename = fullfile(ct_dir, 'coords.csv');
                        if exist(filename, 'file')
                            t_exist = readtableSafe(filename);
                            fprintf('Existing:\n');
                            disp(t_exist);
                            fprintf('New:\n');
                            disp(t);
                            fprintf('An existing file was found here: %s\n', filename);
                            ichoice = menuText('What do you want to do?', {'Add to existing', 'Overwrite and replace (completely)'}, 'multiselect',0);
                            if ichoice == 1
                                vars = intersect(t_exist.Properties.VariableNames, t.Properties.VariableNames);
                                t = [t_exist(:,vars); t(:,vars)];
                            end

                        end
                        writetable(t, filename);
                    end
                    
                end

                % copy ALICE's transform file
                shft = fullfile(zloc, 'ALICE/coregistration/ct_implant_shft.1D');
                aff  = fullfile(zloc, 'ALICE/coregistration/ct_implant_res_shft_al_mat.aff12.1D');
                dest = fullfile(locDirs.ct_1_xfm, 'transform.1D');
                system(sprintf('cat_matvec -ONELINE %s %s > %s', aff, shft, dest));
                assert(exist(fullfile(locDirs.ct_1_xfm, 'transform.1D'),'file') > 0, 'There was a problem copying ALICE transform file to %s\n', dest);
            
                % get leads
                fname = fullfile(ct_dir, 'coords.csv');
                leads = readtableSafe(fname); %SJ - changed from readtable
                
            end % run_xfm
        end

        % -----------------------------------------
        %% ---------- Create Anchors --------------
        % -----------------------------------------
        % Todo: create a matlab GUI to avoid Slicer dependency
        grids = getGrids(subj, rootEEGdir);
        isAnyGrids = ~isempty(grids);
        
        if isAnyGrids
            fprintf('Note to user, the grids are as follows:\n');
            disp(grids);
            if ~exist(fullfile(locDirs.anchors, 'anchors.csv'), 'file') || RERUN_ANCHORS
                anchor_table = anchor(fullfile(locDirs.fs, subj), fullfile(locDirs.mr_pre,'mr_pre.nii'), locDirs.anchors,'slicer_binary',SLICER_BINARY,'postprocess_only',0);
                if ~isempty(anchor_table)
                    writetable(anchor_table, fullfile(locDirs.anchors, 'anchors.csv')); 
                end
            else
                anchor_table = readtableSafe(fullfile(locDirs.anchors, 'anchors.csv')); %SJ - changed from readtable
            end
            
            % utahs
            anchor_tags = anchor_table.tagName;
            anchor_utah_mask = cellfun(@(s) contains(upper(s), 'UTAH'), anchor_tags);
            info_utah_mask = strcmpi(element_info.hardwareType, 'utah');
            if sum(anchor_utah_mask) > 0
                fprintf('Utah arrays detected:\n');
                disp(anchor_tags(anchor_utah_mask));
                if sum(info_utah_mask) ~= numel(unique(anchor_utah_mask))
                    fprintf('Error: The number of utahs in element_info does not match the number you gave in anchors!\n');
                    fprintf('\tYou need to add them to element_info.csv now...\n');
                    unix(sprintf('open %s', fullfile(rootEEGdir, subj, 'docs', 'element_info.csv')));
                    inputYN('Press Enter to continue');
                    
                    utah_table = anchor_table(anchor_utah_mask,:);
                    anchor_table = anchor_table(~anchor_utah_mask,:);
                    fname_utahs = fullfile(locDirs.anchors, 'utah.csv');
                    if ~isempty(utah_table)
                        fprintf('Saving utahs to %s:\n', fname_utahs);
                        disp(utah_table);
                        writetable(utah_table, fname_utahs);
                        writetable(anchor_table, fullfile(locDirs.anchors, 'anchors.csv')); 
                    end
                end
            elseif sum(info_utah_mask) > 0
                warning('Utahs listed in element_info but not marked with Slicer in anchors. Break and rerun or save them yourself separately as utah.csv');
            end
            
            
            
            disp(anchor_table);
            useAnchors = ~isempty(anchor_table) && inputYN('Use these anchors?');
            if ~useAnchors
                anchor_table = [];
            end
        else
            anchor_table = [];
        end

        
        % --------------------------------------------
        %% ---------- Project & save coords ----------
        % --------------------------------------------

        all_coords      = [];
        all_dural_dist  = [];
        all_dural_coord = [];
        all_ct_coords   = [];
        for hemisphere = hems

            fprintf('\n----------------------------------------------------------------------------------------------------\n');
            fprintf('Running hemisphere %s of %s\n', hemisphere, (vector(hems))');
            fprintf('----------------------------------------------------------------------------------------------------\n');
            
            % get the dural and pial surfaces
            dural_filename = fullfile(locDirs.fs_subj,'SUMA',[hemisphere 'h.pial-outer-smoothed.gii']);
            assert(exist(dural_filename, 'file') > 0, 'File not found: %s (please run create_surf.m and suma.m)\n', dural_filename);
            dural_surf = gifti(dural_filename);

            % shouldn't need to do this if it exists already in backup
            if ~USE_ALICE && (run_xfm || RERUN_XFM )
                
                xfm_dir         = locDirs.(sprintf('ct_%d_xfm', ctNum));
                ct_dir          = locDirs.(sprintf('ct_%d', ctNum)); % uigetfile_if_not_exist(default_path, filespec, title, sym, isdir, isBatch)
                ct_coords       = readtableSafe(uigetfile_if_not_exist(fname_default_ct_coord, '*.csv', 'Select CT coordinates table (ct_implant.vox.csv)', 0, 0, BATCH)); %SJ - changed from readtable
                xfm_file        = fullfile(xfm_dir, 'transform.1D');
                if exist(xfm_file, 'file')
                    fprintf('Transform file found; using it: %s\n', xfm_file);
                    % not sure if we should remove parent dir from variable here
                else
                    xfm_file = [];
                end
                xfm_file = register_CT(fullfile(locDirs.mr_pre,'mr_pre.nii'), fullfile(ct_dir, 'ct_implant.nii'), xfm_dir, 'xfm_override', xfm_file);
                if isempty(fileparts(xfm_file))
                    xfm_file = fullfile(xfm_dir, xfm_file);
                end
                
                % quick handle of old csv format (tag + absElmtIdx == chan)
                if ismember('absElmtIdx', ct_coords.Properties.VariableNames) ...
                        && ismember('tagName', ct_coords.Properties.VariableNames) ...
                        && ~ismember('chanName', ct_coords.Properties.VariableNames)
                    ct_coords.chanName = strcat(ct_coords.tagName, arrayfun(@num2str, ct_coords.absElmtIdx, 'UniformOutput',false));
                end

                % Create the afni brik/head files if we don't have it
                ct_nii = fullfile(ct_dir, 'ct_implant.nii');
                ct_afni= fullfile(ct_dir, 'ct_implant+orig.BRIK');
                if exist(ct_nii, 'file') && ~exist(ct_afni, 'file')
                    pwd_cur = pwd;
                    cd(ct_dir)
                    system('3dcopy ct_implant.nii ct_implant'); 
                    cd(pwd_cur);
                end
                
                leads = xfm_leads(ct_coords, xfm_file, fullfile(ct_dir, 'ct_implant+orig'), 'VOXEL', xfm_dir);
                all_ct_coords = cat(1, all_ct_coords, leads);
            end

            projChans   = getLeads(subj, rootEEGdir, 'whichHemi', [hemisphere 'h'], 'hardwareType', {'subdural','micro-subdural'});
            depthChans  = getLeads(subj, rootEEGdir, 'whichHemi', [hemisphere 'h'], 'hardwareType', {'depth'});
            projleads   = leads(ismember(leads.chanName, projChans), :);
            depthleads  = leads(ismember(leads.chanName, depthChans), :);
            info_hem    = element_info(strcmpi(element_info.whichHemi, [hemisphere 'h']), :);
            if isempty(projChans)
                warning('No channels to project');
            end
            
            % Make braindata figure for debugging
            if SHOW_FIGURES
                bd = braindata;
                bd.meta.suma_dir = fullfile(locDirs.fs_subj, 'SUMA');
                bd.loadSubject(subj, rootEEGdir, 'hemi',[hemisphere 'h']);
                bd.plot; bd.setOpacity(0.6); bd.clearPoints;
                % color by tag
                tags = unique(util_split_stringnum(projleads.chanName));
                [~,ndxs] = groupBySubstring(projleads.chanName, tags);
                colors = hsv(numel(tags));
                xyz_mr = projleads{:, {'x','y','z'}};
                for i = 1:numel(ndxs)
                    chan = projleads.chanName(ndxs{i});
                    [~,num] = util_split_stringnum(chan);
                    bd.plotPoint(xyz_mr(ndxs{i},:), 'color', colors(i,:), 'legend', tags{i}, 'label',arrayfun(@(x) {x}, num));
                end
                bd.legend;
                title(sprintf('%s CT %d', subj, ctNum)); 
                set(gcf,'name', sprintf('%s CT %d', subj, ctNum));
                %figfmt;
                pause(5);
            end
            
            if ~isempty(anchor_table)
                if ~ismember('tagName',anchor_table.Properties.VariableNames)
                    anchor_table.tagName = util_split_stringnum(anchor_table.chanName);
                end
                anchor_table_hem = anchor_table(ismember(anchor_table.tagName, info_hem.tagName), :);
                if isempty(anchor_table)
                    warning('Anchor table is empty')
                end
            else
                anchor_table_hem = anchor_table;
            end
            
            
            
            % Run projection spring algorithm
            fname_log = fullfile(locDirs.logs, 'projection');
            fname_fmincon = fullfile(locDirs.logs, 'fmincon');
            [proj_coords, dural_coords, dural_dist] = projection(projleads, dural_surf, anchor_table_hem, info_hem, 'log_file', fname_log, 'fmincon_file',fname_fmincon);
            all_dural_dist = [all_dural_dist; dural_dist];

            % 6/2020 It seems like matlab 2020a has an issue with vertcat for empty tables
            % all_coords = [all_coords; proj_coords; depthleads];
            
            if isempty(all_coords)
                if isempty(proj_coords)
                    all_coords = depthleads;
                elseif isempty(depthleads)
                    all_coords = proj_coords;
                else
                    all_coords = [proj_coords; depthleads];
                end
            else
                if isempty(proj_coords)
                    all_coords = [all_coords; depthleads];
                elseif isempty(depthleads)
                    all_coords = [all_coords; proj_coords];
                else
                    all_coords = [all_coords; proj_coords; depthleads];
                end
            end
            
            %all_dural_coord = [all_dural_coord; dural_coords; depthleads];
            if isempty(all_dural_coord)
                if isempty(dural_coords)
                    all_dural_coord = depthleads;
                elseif isempty(depthleads)
                    all_dural_coord = dural_coords;
                else
                    all_dural_coord = [dural_coords; depthleads];
                end
            else
                if isempty(dural_coords)
                    all_dural_coord = [all_dural_coord; depthleads];
                elseif isempty(depthleads)
                    all_dural_coord = [all_dural_coord; dural_coords];
                else
                    all_dural_coord = [all_dural_coords; dural_coords; depthleads];
                end
            end
            
            % debugging figure
            if SHOW_FIGURES && ~isempty(proj_coords)
                bd = braindata;
                bd.meta.suma_dir = fullfile(locDirs.fs_subj, 'SUMA');
                bd.loadSubject(subj, rootEEGdir, 'hemi',[hemisphere 'h']);
                bd.plot; bd.setOpacity(0.6); bd.clearPoints;
                % color by tag
                tags = unique(util_split_stringnum(proj_coords.chanName));
                [~,ndxs] = groupBySubstring(proj_coords.chanName, tags);
                colors = hsv(numel(tags));
                xyz = proj_coords{:, {'x','y','z'}};
                for i = 1:numel(ndxs)
                    chan = proj_coords.chanName(ndxs{i});
                    [~,num] = util_split_stringnum(chan);
                    bd.plotPoint(xyz(ndxs{i},:), 'color', colors(i,:), 'legend', tags{i}, 'label',arrayfun(@(x) {x}, num));
                end
                bd.legend;
                title(sprintf('%s Projected (CT %d)', subj, ctNum));
                set(gcf,'name', sprintf('%s Projected (CT %d)', subj, ctNum));
                %figfmt;
                pause(5);
            end
            
        end % hemisphere
        
        fname_coords = fullfile(locDirs.mr_pre, sprintf('coords_%d.csv', ctNum));
        fname_snap   = fullfile(locDirs.mr_pre, sprintf('coords_snap_%d.csv', ctNum));
        fname_ddist  = fullfile(locDirs.mr_pre, sprintf('dural_dist_%d.csv', ctNum));
        %fname_ctcoords = fullfile(locDirs.mr_pre, sprintf('ct_coords_%d.csv', ctNum));
        writetable(all_coords, fname_coords);
        writetable(all_dural_coord, fname_snap);
        writetable(all_dural_dist, fname_ddist);
        %writetable(all_ct_coords, fname_ctcoords);
    end % CTNum
    
    % End of projection-per-CT coordinate section
    
    
    % --------------------------------------------
    %% Compile coords_1-2.csv into leads.csv
    % --------------------------------------------
    fname = fullfile(locDirs.mr_pre, 'coords_1.csv');
    all_coords = readtableSafe(fname); %SJ - changed from readtable
    if numCTs == 2
        [tagRemapTable, chan_summary, tag_summary] = reconcileMultiCT(rootEEGdir, subj, locDirs, element_info, BATCH);
        fname_chan_summary = fullfile(locDirs.mr_pre, 'shift_summary_chan.csv');
        fname_tag_summary  = fullfile(locDirs.mr_pre, 'shift_summary_tag.csv');
        fname_rename_summary= fullfile(locDirs.mr_pre, 'shift_rename.csv');
        writetable(chan_summary, fname_chan_summary);
        writetable(tag_summary, fname_tag_summary);
        writetable(tagRemapTable, fname_rename_summary);
        
        if SHOW_FIGURES
            bp = brainplotter;
            figure;
            bd = braindata2(subj, rootEEGdir);
            for h=hems
                bp.loadSurface(bd, ['pial ' char(h) 'h']);
            end
            bp.plot(fieldnames(bp.surfaces));
            bp.plotPoint(bd.zloc.proj_1_xyz, 'color',[0 1 0], 'legend','CT_1');
            bp.plotPoint(bd.zloc.proj_2_xyz, 'color',[0 0 1], 'legend','CT_2');
            bp.legend('orientation','horizontal');
            
        end
        
        
        fname2 = fullfile(locDirs.mr_pre, 'coords_2.csv');
        t1 = readtableSafe(fname); %SJ - changed from readtable
        t2 = readtableSafe(fname2); %SJ - changed from readtable
        vars = intersect(t1.Properties.VariableNames, t2.Properties.VariableNames);
        t1 = t1(:,vars);
        t2 = t2(:,vars);
        
        [tagName2, chanNum2] = util_split_stringnum(t2.chanName);
        newName = ismember(tagName2, tagRemapTable.old_tag);
        update  = ~newName;
        newName = find(newName);
        update  = find(update);
        
        % update (or add new)
        for i = 1:numel(update)
            % for each row in t2...
            chan = t2.chanName(update(i));
            if ismember(chan, t1.chanName)
                % replace the CT1 row
                t1(strcmpi(t1.chanName, chan), :) = t2(strcmpi(t2.chanName, chan), :);
            else
                % add completely new channel as new row
                row = t2(update(i), :);
                row.chanName = util_combine_strnum(tagName2(update(i)), chanNum2(update(i)));
               % fprintf('%s - %d\n', char(row.chanName), i);
                t1 = [t1; row];
            end
        end
        
        % add old channels which have shifted as if they were new channels
        for i = 1:numel(newName)
            row = t2(newName(i), :);
            temp = tagName2(newName(i));
            newTag = tagRemapTable(strcmpi(tagRemapTable.old_tag, temp), :).rename_tag;
            row.chanName = util_combine_strnum(newTag, chanNum2(newName(i)));
            t1 = [t1; row];
        end
        
        all_coords = t1(:,{'chanName','x','y','z'});
        all_coords{:,{'x','y','z'}} = round(all_coords{:,{'x','y','z'}}, 2);
        fname = fullfile(locDirs.mr_pre, 'leads.csv');
        assert(numel(all_coords.chanName) == numel(unique(all_coords.chanName)), 'Channel names are not unique after ammendment');
        writetable(all_coords, fname);
        
    else
        fname = fullfile(locDirs.mr_pre, 'leads.csv');
        all_coords{:,{'x','y','z'}} = round(all_coords{:,{'x','y','z'}}, 2);
        all_coords = all_coords(:,{'chanName','x','y','z'});
        writetable(all_coords, fname);
    end
    
    % Compute and save the bipolar pseudo-electrodes based on euclidean midpoint
    createBipolar(subj, rootEEGdir);
    
    
    % -----------------------------------------
    %% ---------- Create ROIs ----------

    localizer_rois(subj, rootEEGdir, RERUN_ROI);
    localizer_atlas(subj, rootEEGdir, RERUN_ATLAS);
    
    % --------------------------------------------
    %% Write files from zloc into tal
    % --------------------------------------------
    copyfile(fullfile(locDirs.mr_pre, 'leads.csv'), tal);
    fprintf('leads.csv written to tal\n');
    
    checkTal(subj, rootEEGdir);
end





function empty = is_empty_dir(filepath)
    s = dir(filepath);
    charname = char(s.name);
    count = sum(charname(:,1) ~= '.');
    empty = count == 0;
end

function copyOrigs(rootEEGdir, locDirs, subj, TAL_ORIG, rootEEGdir_ORIG, USE_ALICE, numCTs, TAL_ORIG_SUB)
    fprintf('Copying files to the rootEEGdir directory, this may take a minute...');
    
    if  exist('TAL_ORIG_SUB','var')
        dsubj_talOrig = fullfile(TAL_ORIG, subj, TAL_ORIG_SUB);
    else
        dsubj_talOrig = fullfile(TAL_ORIG, subj);
    end
    
     % docs
%     if ~exist(fullfile(locDirs.docs),'dir') || ~exist(fullfile(locDirs.docs, 'element_info.csv'), 'file')
%         mkdir(locDirs.docs);
%         copyfile(fullfile(rootEEGdir_ORIG,subj,'docs'), locDirs.docs); 
%     end

    % CT coords
    if is_empty_dir(locDirs.ct_1) 
        orig = fullfile(dsubj_talOrig,'/intermediates/locs_0_curry/coords.ct_implant.vox.csv');
        if exist(orig, 'file')
            copyfile(orig, locDirs.ct_1);
        end
    end
    if numCTs == 2 && is_empty_dir(locDirs.ct_2)
        orig = fullfile(dsubj_talOrig,'/intermediates/locs_0_curry/coords.ct_implant_2.vox.csv');
        if exist(orig, 'file')
            copyfile(orig, locDirs.ct_2);
        end
    end

    % freesurfer
    if is_empty_dir(locDirs.fs_subj) && exist(fullfile(dsubj_talOrig,'/intermediates/images_2_fsSumaStd'), 'file')
        copyfile(fullfile(dsubj_talOrig,'/intermediates/images_2_fsSumaStd'), locDirs.fs);
    end

    % MR dicom / nii
    if is_empty_dir(locDirs.mr_pre_dicom)
        orig = fullfile(dsubj_talOrig,'/intermediates/images_0_dicom/mr_pre');
        if exist(orig, 'file')
            copyfile(orig, locDirs.mr_pre_dicom);
        else
            % try .nii
            nii = fullfile(locDirs.mr_pre, 'mr_pre.nii');
            orig = fullfile(dsubj_talOrig,'/intermediates/images_1_anon/mr_pre/mr_pre.nii');
            if ~exist(nii, 'file')
                copyfile(orig, nii);
            end
        end
    end

    % CT dicom /nii
    if is_empty_dir(locDirs.ct_1_dicom)
        orig = fullfile(dsubj_talOrig,'/intermediates/images_0_dicom/ct_implant');
        if exist(orig, 'file')
            copyfile(orig, locDirs.ct_1_dicom);
        else
            % try .nii
            nii = fullfile(locDirs.ct_1, 'ct_implant.nii');
            orig = fullfile(dsubj_talOrig,'/intermediates/images_1_anon/ct_implant/ct_implant.nii');
            if ~exist(nii, 'file')
                copyfile(orig, nii);
            end
        end
    end
    % CT dicom /nii 2
    if numCTs == 2 && is_empty_dir(locDirs.ct_2_dicom)
        orig = fullfile(dsubj_talOrig,'/intermediates/images_0_dicom/ct_implant_2');
        if exist(orig, 'file')
            copyfile(orig, locDirs.ct_2_dicom);
        else
            % try .nii
            nii = fullfile(locDirs.ct_2, 'ct_implant_2.nii');
            orig = fullfile(dsubj_talOrig,'/intermediates/images_1_anon/ct_implant_2/ct_implant_2.nii');
            if ~exist(nii, 'file')
                copyfile(orig, nii);
            end
        end
    end

    % anchors
    if is_empty_dir(locDirs.anchors) && ~isempty(getGrids(subj, rootEEGdir))
        anchor_orig = fullfile(dsubj_talOrig,'/intermediates/locs_0_slicer/anchors.csv');
        if exist(anchor_orig, 'file')
            copyfile(anchor_orig, locDirs.anchors);
        else
            fcsv_orig = fullfile(dsubj_talOrig,'/intermediates/locs_0_slicer/anchors.fcsv');
            if exist(fcsv_orig, 'file')
                slicer_fcsv2csv(fcsv_orig);
                copyfile(anchor_orig, locDirs.anchors);
            end
        end
    end

    % Transform
    if ~USE_ALICE && is_empty_dir(locDirs.ct_1_xfm)
        dreg = fullfile(dsubj_talOrig,'/intermediates/locs_1_registered/ict2bmrAF');
        dxfm = fullfile(dsubj_talOrig,'/intermediates/images_3_xfm/ict2bmrAF');
        assert(exist(dreg,'dir') > 0);
        assert(exist(dxfm,'dir') > 0);
        ct2marker = '_2'; % marker for 2nd CT

        % registered/ict2bmrAF
        files = lsCell(dreg)';
        for i = 1:numel(files)
            name = files{i};
            if contains(name, 'csv'), continue; end
            if strfound(name, ct2marker) && numCTs == 2
                % ct_2
                if strfound(name, 'coords.mr_pre_2.ict2bmrAF.gii') % we renamed this simply 'elec'
                    copyfile(fullfile(dreg, name), fullfile(locDirs.ct_2_xfm,'elec.gii'));
                elseif ~strfound(name, 'elec.gii')
                    copyfile(fullfile(dreg, name), fullfile(locDirs.ct_2_xfm, strrep(name,ct2marker,'')));
                end
            elseif ~strfound(name, ct2marker)
                % ct_1
                if strfound(name, 'coords.mr_pre.ict2bmrAF.gii') % we renamed this simply 'elec'
                    copyfile(fullfile(dreg, name), fullfile(locDirs.ct_1_xfm,'elec.gii'));
                elseif ~strfound(name, 'elec.gii')
                    copyfile(fullfile(dreg, name), fullfile(locDirs.ct_1_xfm, name));
                end
            end
        end

        % coord files
        file = fullfile(dreg, '../coords.mr_pre.ict2bmrAF.csv');
        if exist(file,'file'), copyfile(file, fullfile(locDirs.ct_1_xfm, 'coords.csv')); end
        file = fullfile(dreg, '../coords.mr_pre_2.ict2bmrAF.csv');
        if exist(file,'file'), copyfile(file, fullfile(locDirs.ct_2_xfm, 'coords.csv')); end

        files = lsCell(dxfm)';
        for i = 1:numel(files)
            name = files{i};
            if strfound(name, ct2marker) && numCTs == 2
                % ct2
                if strfound(name, 'img_check_color')
                    copyfile(fullfile(dxfm,name), locDirs.ct_2_color);
                elseif strfound(name, 'img_check_edge')
                    copyfile(fullfile(dxfm,name), locDirs.ct_2_edge);
                else
                    copyfile(fullfile(dxfm,name), fullfile(locDirs.ct_2_xfm, strrep(name,ct2marker,'')));
                end
            elseif ~strfound(name, ct2marker)
                % ct1
                if strfound(name, 'img_check_color')
                    copyfile(fullfile(dxfm,name), locDirs.ct_1_color);
                elseif strfound(name, 'img_check_edge')
                    copyfile(fullfile(dxfm,name), locDirs.ct_1_edge);
                else
                    copyfile(fullfile(dxfm,name), fullfile(locDirs.ct_1_xfm, name));
                end
            end
        end
    end
    fprintf('done!\n');
end

function ctmr = run_alice(zloc, subj, locDirs, ctNum, hemisphere)
    rootEEGdirEEG = fullfile(zloc, '../../..'); %rootEEGdir/subj/tal/zloc

    f = findall(groot, 'name', 'ALICE: Also Developed by a Cast of Thousands');
    close(f);
    
    ctmr = ctmrGUI();
    dalice = fullfile(zloc, 'ALICE');
    
    % Load/create the ALICE directory
    cd(zloc);
    if ~exist(dalice,'dir')
        ctmr.CreateDirectory();
    else
        ctmr.LocateDirectory(dalice);
    end
    cd(zloc);
    
    % set name
    set(ctmr.controls.edtSbjName, 'String',subj);
    ctmr.settings.subject = subj;
    
    % set MRI
    if ~exist(fullfile(dalice,'MRI'),'file'), mkdir(fullfile(dalice,'MRI')); end
    % Recently cases have not been aligning from their MR, so use the freesurfer-->suma MR volume
    % ^ Update: actually, the MR just needs to be de-obliqued!
    mr = fullfile(locDirs.fs_subj, 'SUMA', sprintf('%s_SurfVol.nii',subj));
    fprintf('Loading MR: %s\n', mr);
    ctmr.btnOpenMRI([],[],mr);

    % set FS
    if ~exist(fullfile(dalice,'FreeSurfer'),'file'), mkdir(fullfile(dalice,'FreeSurfer')); end
    unix(sprintf('cp %s %s', fullfile(locDirs.fs_subj, 'mri/ribbon.mgz'), fullfile(locDirs.fs_subj,'SUMA')));
    fs      = fullfile(locDirs.fs_subj, 'SUMA', 'ribbon.nii');
    fs_mgz  = fullfile(locDirs.fs_subj, 'SUMA', 'ribbon.mgz');
    if ~exist(fs, 'file') && exist(fs_mgz, 'file')
        pwd_cur = pwd;
        cd(fullfile(locDirs.fs_subj, 'SUMA'));
        system('mri_convert ribbon.mgz ribbon.nii');
        cd(pwd_cur);
    else
        fprintf('Could not find a 2-hemisphere ribbon (%s) using single hemisphere\n', fs_mgz);
        fs = fullfile(locDirs.fs_subj, 'SUMA', sprintf('%sh.ribbon.nii',hemisphere));
        assert(exist(fs, 'file') > 0, 'No ribbons found in SUMA. Check for SUMA/Freesurfer files/issues');
    end
    fprintf('Loading ribbon: %s\n', fs);
    ctmr.btnOpenFS([],[],fs);
    
    % set CT
    if ~exist(fullfile(dalice,'CT'),'file'), mkdir(fullfile(dalice,'CT')); end
    if ctNum == 1, dct = locDirs.ct_1; else, dct = locDirs.ct_2; end
    ct = fullfile(dct, 'ct_implant.nii');
    fprintf('Loading CT: %s\n', ct);
    ctmr.btnOpenCT1([],[],ct);

    fprintf('IF you want to coregister BEFORE clustering, select the coregistered CT before STEP 2\n');
    
    root = fullfile(zloc,'../../..');
    % set hemisphere
    if hemisphere == 'r'
        ctmr.controls.radiobtn3.Value = 0;
        ctmr.controls.radiobtn4.Value = 1;
        ctmr.radiobtnSelectionHemisphere([], struct('NewValue', struct('String', 'Right')));
        chans = getLeads(subj, root, 'whichHemi', 'rh', 'implant',ctNum);
        
    else % left
        ctmr.controls.radiobtn3.Value = 1;
        ctmr.controls.radiobtn4.Value = 0;
        ctmr.radiobtnSelectionHemisphere([], struct('NewValue', struct('String', 'Left')));
        chans = getLeads(subj, root, 'whichHemi', 'lh', 'implant',ctNum);
    end
    
    if isempty(chans)
        error('getLeads returned no chans for %s %s %s implant %d', subj, root, hemisphere, ctNum);
    end
    
    if numel(unique(chans)) ~= numel(chans)
        error('channels are not all unique');
    end
    
    ctmr.controls.radiobtn3.Enable = 'off';
    ctmr.controls.radiobtn4.Enable = 'off';
    
    ctmr.settings.chans = chans;
  
    choice = [];
    saved = 0;
    options = {
        'Keyboard' % 1
        'Clear previous clusters' % 2
        'Electrode number/name table (during selection)'
        'Relaunch GUI' % 4
        'Plot / Save Coordinates (after selection)'
        'View previously extracted cluster settings' % 6
        'Continue / Next hemisphere'
        'Copy Updated Afni scripts from Source' % 8
        };
    while isempty(choice) || choice ~= 7 || ~saved
        if choice == 7
            if inputYN('Careful--you have not yet saved the electrodes yet, continue anyway?')
                fprintf('Setting saved flag to true. Menu will display again, this time allowing your continue\n')
                saved = 1;
            end
        end
        
        choice = menuText(['Use GUI to cluster/select electrodes: ' hemisphere], options{:}, 'multiSelect',0);
        switch choice
            case 9
                % TODO: mask the skull. reduction in ROIs might improve speed
                
                cd(fileparts(ctmr.settings.MRI));
                fprintf('Calculating the skull strip. This may take a few minutes....');
                command = sprintf([ ...
                            '3dSkullStrip -input %s_SurfVol.nii -prefix mr.ns.nii; ' ...
                            '3dmask_tool -input mr.ns.nii -dilate_input 5; ' ...                
                            '3dcalc -a combined_mask+orig -b ../coregistration/ct_implant_res_shft_al+orig -expr "a*b" -overwrite;' ...
                            '3dcopy calc+orig ct_implant_xfmd_ns.nii' ], ...
                            subj);

                system(command);
                cd(ctmr.settings.currdir);
                fprintf('done!\n');
                fprintf('Now TODO: use ct_implant_xfmd_ns.nii instead of transformed CT\n');
               
            case 8
                % Copy updated scripts from root ALICE3 directory
                ctmr.settings.scriptspath  = [fileparts( mfilename('fullpath') ) '/'];
                
                %locate afni scripts
                %dalice
                obj.settings.scriptspath  = [fileparts( mfilename('fullpath') ) '/'];
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/alignCTtoT1_shft_res.csh'], [dalice '/coregistration/']);
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/3dclustering.csh'], [dalice '/3Dclustering/']);
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/select_electrode.csh'], [dalice '/3Dclustering/']);
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/open_afni_suma.csh'], [dalice '/3Dclustering/']);
                copyfileSafe([obj.settings.scriptspath 'AFNI_scripts' '/indexify_electrodes.csh'], [dalice '/3Dclustering/']);
                        
            case 1
                % You may stay at this keyboard as long as you'd like. The GUI runs in a separate thread
                disp('At keyboard. Use run button to continue.');
                keyboard
            case 2 % Clear cluster files to allow reclustering
                system(sprintf('mv %s ~/.Trash', fullfile(dalice,'3Dclustering/3dclusters*')));
            case 3
                if isempty(ctmr.subGUI)
                    disp('Option can only be used while selecting electrodes')
                else
                    fprintf('%s Hemisphere (x=selected):\n', ctmr.settings.Hemisphere);
                    for i = 1:numel(chans)
                        if ismember(i,ctmr.subGUI.selected_list), sel='x'; else sel=''; end
                        fprintf('%1s %2d - %s\n', sel , i, chans{i});
                    end
                end
                
            case 4
                ctmr = run_alice(zloc, subj, locDirs, ctNum, hemisphere);
                return;
            
            case 6
                fprintf('\nExisting extracted cluster files:')
                disp(char(lsCell(fullfile(dalice,'3Dclustering/3dclusters*.gii'))));
                fprintf('(filename key: r=volume, is=spacing, thr=intensity)\n');
            
            case 5 % Plot/Save
            
                % transforma coords to MR coordinate system
                % extract CM using AFNI-SUMA plug-in.
                CM = importdata(fullfile(dalice,'/3Dclustering/electrode_CM.1D'));
                CM.data = flip(CM.data); % --> most recent first
                
                %remove repeated electrodes by taking a unique starting with the most recent
                [~, index] = unique(CM.data(:,4), 'stable');
                n = min(numel(chans), numel(index));
                if n < numel(chans)
                    warning('There are %d channels but only %d marked in CT', numel(chans), n);
                end
                
                index = index(1:n);

                % put in ascending order
                CM = CM.data(index,1:4);
                [~,order] = sort(CM(:,4));

                elecCoord  = CM(:,1:3);
                elecNum    = CM(:,4);       % corresponds to chans

                %check for empty rows and put NANs
                xyz_ct = CM(:,1:3); % create empty array 

                USE_LEGACY_XFM          = 0;
                CT_IN_MR_SPACE_ALREADY  = 0;
                
                
                if USE_LEGACY_XFM == 1
                    % write to file
                    
                    these_chans = [];
                    for i = 1:numel(elecNum)
                        these_chans{i} = chans{elecNum(i)};
                        fd = fopen(fullfile(locDirs.ct_1_xfm, 'mm.txt'),'a');
                        fprintf(fd, '%f %f %f\n', xyz_ct(i,1), xyz_ct(i,2), xyz_ct(i,3));
                        fclose(fd);
                    end
                    ct_coords = table(these_chans',xyz_ct(:,1),xyz_ct(:,2),xyz_ct(:,3),'VariableNames',{'chanName','x','y','z'});
                    
                    affine_xfm = fullfile(locDirs.ct_1_xfm,'CT_highresRAI_res_shft_al_mat.aff12.1D');
                    shift_xfm = fullfile(locDirs.ct_1_xfm,'CT_highresRAI_shft.1D');
                    fname_xfm = fullfile(locDirs.ct_1_xfm, 'transform.1D');
                    command = sprintf('cat_matvec -ONELINE %s %s > %s', affine_xfm, shift_xfm, fname_xfm);
                    unix(command);
                    mr = xfm_leads(ct_coords, fname_xfm, fullfile(locDirs.ct_1, 'ct_implant+orig'), 'RAI');
                
                elseif USE_LEGACY_XFM == 2
                    ct_coords = table(these_chans',xyz_ct(:,1),xyz_ct(:,2),xyz_ct(:,3),'VariableNames',{'chanName','x','y','z'});
                    fname_xfm = fullfile(locDirs.ct_1_xfm,'transform.1D');
                    assert(exist(fname_xfm, 'file') > 0, 'Not found: %s', fname_xfm);
                    mr = xfm_leads(ct_coords, fname_xfm, fullfile(locDirs.ct_1, 'ct_implant+orig'), 'RAI');
                end
                    
                
                if CT_IN_MR_SPACE_ALREADY
                    xyz_mr = xyz_ct;
                    xyz_mr = [-xyz_mr(:,1:2) xyz_mr(:,3)]; 
                else
                    
                    % Get tranform matrix files
                    coreg_files = lsCell(fullfile(dalice,'coregistration'));
                    fname_aff   = fullfile(dalice, 'coregistration', coreg_files(contains(coreg_files, 'aff12')));
                    if numel(fname_aff) > 1
                        fprintf('You usually want the ct_implant_res_shft_al_mat one\n');
                        [~,fname_aff] = menuText('More than 1 affine file, choose:', fname_aff, 'multiselect',0);
                    else
                        fname_aff = char(fname_aff);
                    end
                    
                    fname_shft  = fullfile(dalice, 'coregistration', coreg_files(contains(coreg_files, 'shft.1D')));
                    if numel(fname_shft) > 1
                        [~,fname_shft] = menuText('More than 1 shift file, choose:', fname_shft, 'multiselect',0);
                    elseif isempty(fname_shft)
                        error('No shft.1D file found in coregistration folder');
                    else
                        fname_shft = char(fname_shft);
                    end
                    
                    
                    % apply rigid xfm, then affine xfm, then go from LPI to RAI by negating x/y coord
                    Taffin = [reshape(dlmread(fname_aff), [4 3])'; [0 0 0 1]];
                    Tshift = [reshape(dlmread(fname_shft), [4 3])'; [0 0 0 1]];
                    Tfull = Tshift * Taffin;
                    X = cat(2,xyz_ct, ones(size(xyz_ct,1),1));
                    xyz_mr = (inv(Tfull) * X')';
                    xyz_mr = [-xyz_mr(:,1:2) xyz_mr(:,3)]; 
                end
                
                
                bd = braindata;
                bd.meta.suma_dir = fullfile(locDirs.fs_subj, 'SUMA');
                bd.loadSubject(subj, rootEEGdirEEG, 'hemi',[hemisphere 'h']);
                bd.plot; bd.setOpacity(0.6); bd.clearPoints;
                % color by tag
                these_chans = chans(elecNum);
                tags = unique(util_split_stringnum(these_chans));
                [~,ndxs] = groupBySubstring(these_chans, tags);
                colors = hsv(numel(tags));
                for i = 1:numel(ndxs)
                    chan = these_chans(ndxs{i});
                    [~,num] = util_split_stringnum(chan);
                    bd.plotPoint(xyz_mr(ndxs{i},:), 'color', colors(i,:), 'legend', tags{i}, 'label',arrayfun(@(x) {x}, num));
                end
                bd.legend;
                %figfmt;

                % display/save to a hemisphere file
                for i = 1:numel(elecNum)
                    chan = chans{elecNum(i)};
                    fprintf('%f %f %f %s\n', xyz_mr(i,1), xyz_mr(i,2), xyz_mr(i,3), chan);
                end

                ct_xfm_dir = locDirs.(sprintf('ct_%d_xfm', ctNum));
                
                if inputYN('Save these coordinates?')
                    % Note we append (not overwrite)
                    fd = fopen(fullfile(ct_xfm_dir, sprintf('coords_%s.txt', hemisphere)), 'a');
                    for i = 1:numel(elecNum)
                        chan = chans{elecNum(i)};
                        fprintf(fd, '%f %f %f %s\n', xyz_mr(i,1), xyz_mr(i,2), xyz_mr(i,3), chan);
                    end
                    fclose(fd);
                    saved = 1;
                else
                    disp('If things look very wrong, one possible reason is that you used a different max intensity');
                    disp(' between different selections. This is not allowed and things will be wrong on your plot.');
                end
                
            case 7
                % Loop again
            otherwise
                disp('Invalid selection');
        end
                
    end
    
    
end

function [tagRemapTable, chan_summary, tag_summary] = reconcileMultiCT(rootEEGdir, subj, locDirs, element_info, isBatch)
    % This function was recycled from old code and depends on a few old functions
    %
    % Looks at the xyz coordinates based on different CT's. If euclidean
    % difference is over the move_threshold, consider the moved electrode a new
    % electrode. Otherwise, just use the first location
    % Actually, the way we do it is if any electrode on a peice of hardware
    % moves, make new electrodes for all the new hardware
    %
    % Note: it needs to be OK if the two location files are MISSING some leads
    % from the jacksheet!
   
    MOVE_THRESHOLD_AVG  = 2.5; % 2.5 mm (average for hardware)
    MOVE_THRESHOLD_ANY  = 4.5; % 4.0 mm (any electrode on a strip/grid)
    fnameInfoBak    = fullfile(rootEEGdir, subj, 'docs/element_info.csv.bak');
    fnameInfo       = fullfile(rootEEGdir, subj, 'docs/element_info.csv');
    fnameJack       = fullfile(rootEEGdir, subj, 'docs/jacksheetMaster.csv');
    fnameJackBak    = fullfile(rootEEGdir, subj, 'docs/jacksheetMaster.csv.bak');
   
    if nargin < 5, isBatch = 0; end
    
    fname1 = fullfile(locDirs.mr_pre, 'coords_1.csv');
    fname2 = fullfile(locDirs.mr_pre, 'coords_2.csv');
    t1 = readtableSafe(fname1);
    t2 = readtableSafe(fname2);
    oldInfo = element_info;
        
    [tagNames, absElmtIdx] = util_split_stringnum(t1.chanName);
    t1.tagName = tagNames;
    t1.absElmtIdx = absElmtIdx;
    
    [tagNames, absElmtIdx] = util_split_stringnum(t2.chanName);
    t2.tagName = tagNames;
    t2.absElmtIdx = absElmtIdx;
    
    % correcting case of one table having too many columns
    colNames = intersect(t1.Properties.VariableNames, t2.Properties.VariableNames);
    t1 = t1(:, colNames);
    t2 = t2(:, colNames);
    
    % Ensure these are sorted correctly
    tags = intersect(element_info.tagName, union(util_split_stringnum(t1.chanName), util_split_stringnum(t2.chanName)), 'stable');
    [~,presortOrder] = sort(t1.absElmtIdx);
    t1 = t1(presortOrder, :);
    order = sortBySubstring(t1.tagName, tags);
    t1 = t1(order, :);
    [~,presortOrder] = sort(t2.absElmtIdx);
    t2 = t2(presortOrder, :);
    order = sortBySubstring(t2.tagName, tags);
    t2 = t2(order, :);
    
    % Run the combiner scripts
    [comTable, allTable] = ip_compare_lead_xyz(t1, t2); % find differences 
    [~, tagRemapTable, tag_summary] = ip_combiner_v2(comTable, allTable, MOVE_THRESHOLD_ANY, MOVE_THRESHOLD_AVG);
    chan_summary = comTable(:, {'chanName','d_xyz'});
    chan_summary.d_xyz = round(chan_summary.d_xyz, 2);
    tag_summary{:,2:end} = round(tag_summary{:,2:end}, 2);

    disp(chan_summary);
    disp(tag_summary);
    disp(tagRemapTable);
    
    % Are there any new tags that aren't in element_info?
    tagRemapTableSub = tagRemapTable(ismember(tagRemapTable.rename_tag, setdiff(tagRemapTable.rename_tag, oldInfo.tagName)), :);
    if isempty(tagRemapTableSub)
        fprintf('Either there is no need to rename element_info or things have been renamed there already\n');
    else
        %-- update element_info--%
        
        while ~isBatch && ~inputYN('Investigate the above. Okay to proceed and save? ')
            keyboard;
        end
        
        while ~isBatch && ~inputYN('Have you updated the dateOut field? (N - open/edit element_info, Y - continue)')
            system(sprintf('open %s', fnameInfo));
            keyboard;
            oldInfo = getElementInfo(subj, rootEEGdir);
        end

        if isBatch || inputYN(sprintf('Should I update the element_info and then createMasterJack for %s?', subj))
            % update element_info and jacksheet master
            newInfo = oldInfo;
            for i = 1 : height(tagRemapTableSub)
                newTag = tagRemapTableSub{i, 'rename_tag'};
                oldTag = tagRemapTableSub{i, 'old_tag'};
                oldRowNdx = find(strcmpi(newInfo.tagName, oldTag));
                oldRow = newInfo(oldRowNdx, :);
                newRow = oldRow;
                newRow.tagName = newTag;
                newRow.isInterictal = {'[]'};
                newRow.isResected = {'[]'};
                newRow.isIctal = {'[]'};
                newRow.dateIn = oldRow.dateOut;
                newRow.dateOut = {''};
                newRow.notes = {'Shifted element (detected by localizer)'};

                newInfo = [newInfo; newRow];
                newInfo = table_swap(newInfo, height(newInfo), 1 + oldRowNdx);
            end

            % now we have to move the CLIN and SYNC channels to the end....
            ndxs = find(strcmpi(newInfo.chanType, 'CLIN'));
            ndxs = flip(ndxs);
            for i = 1:numel(ndxs)
                % bubble down
                for j = ndxs(i) : height(newInfo)-i
                    newInfo = table_swap(newInfo, j, j+1);
                end
            end
            
            ndxs = find(strcmpi(newInfo.chanType, 'SYNC'));
            ndxs = flip(ndxs);
            for i = 1:numel(ndxs)
                % bubble down
                for j = ndxs(i) : height(newInfo)-i
                    newInfo = table_swap(newInfo, j, j+1);
                end
            end
            
            disp(newInfo);
            if ~isBatch && ~inputYN('Does that look good? Type "no" to make edits yourself')
                system(sprintf('open %s', fnameInfo));
                keyboard;
            else
                % bakckup and write
                if ~exist(fnameInfoBak, 'file')
                    movefile(fnameInfo, fnameInfoBak);
                end
                writetable(newInfo, fnameInfo);
                
                if ~exist(fnameJackBak, 'file')
                    movefile(fnameJack, fnameJackBak);
                end
                createMasterJack(subj, rootEEGdir);
            end
            
           
        end
    end    
end % reconcileMultiCT

function newt = table_swap(t, ri, rj)
% swap ri'th row with rj'th row
    if ri==rj
        newt = t;
        return; 
    end
    
    temp = sort([ri,rj]);
    i = temp(1);
    j = temp(2);

    irow = t(i,:);
    jrow = t(j,:);
    newt = [t(1:i-1,:); jrow; t(i+1:j-1, :); irow; t(j+1:end, :)];
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