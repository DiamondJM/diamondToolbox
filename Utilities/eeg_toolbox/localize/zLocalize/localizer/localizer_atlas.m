function localizer_atlas(subj, rootEEGdir, overwrite)
    % localizer_atlas(subj, root, [overwrite=0]) reformats post-projection native afni/freesurfer atlas files into a more useful format
    %
    % Monopolar and bipolar atlas files will be placed in the appropriate locations:
    %   [subj]/tal/atlas/[mono|bi]polar.mat
    %
    % Output file format details:
    %   - Define electrode "e",  atlas "a", atlas-label "l_a", and electrode-distance-to-label, "d_{e,l_a}", 
    %       and consider all possible unique tuples, (e, a, l_a, d_{e,l_a}). The table row are defined by these tuples.
    %
    %   - This table contains a row for all unique tuples that satisfy:
    %       - d_{e,l_a} <= MAX_DIST
    %       - For any particular (e*,a*) pair, there is at most one row (e*, a*,l_a*, d{e*,l_a*}),
    %           where (l_a*, d{e*,l_a*}) is defined by l_a* = argmin_{l_a} d{e*, l_a}.
    %           (i.e. We find the nearest label for each atlas)
    %
    %   - Labels containing the word "unknown" are not included if there is an alternative label within MAX_DIST   
    %
    %   - Note_2: AFNI has some annoying behavior. For all cortical parcellations (e.g. aparc, and aparc.a2009s), it
    %           creates a segmentation-parcellation "combo", defined by the parcellation overlayed onto the segmentation
    %           and named accordingly (e.g. aparc+aseg, aparc.a2009s+aseg) (note parcellation labels take precedence).
    %
    %     * Since we are often interested in ONLY the parcellation WITHOUT any segmentation labels, we create our own
    %       parc-only atlases (e.g. aparc, aparc.a2009s) defined by the labels in aparc+aseg without the labels in aseg.
    %   
    %     Ex, NIH059:
    %       - (TT1, aseg+aparc, Left-Cerebellum-Cortex *, 0 mm)
    %       - (TT1, aseg+aparc, ctx-lh-fusiform *, 1mm)             >> (TT1, aparc, ctx-lh-fusiform *, 1mm)
    %
    %       * Left-Cerebellum-Cortex is from the segmentation atlas, 
    %         while ctx-lh-fusiform is from the parcellation atlas
    %       
    %       >> We want (TT1, aparc, ctx-lh-fusiform *, 1mm) to be in our table too, even though Left-Cerebellum-Cortex
    %          is closer in the aparc+aseg "combo" atlas.
    %
    % Resources:
    %     Freesurfer atlas outputs: 
    %         https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT
    %
    %     Destrieux parc table:
    %         "Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature"
    %         https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2937159/
    %         $FREESURFER_HOME/Simple_surface_labels2009.txt 
    %
    % REVISION_HISTORY
    %   ??/?? MST - Created
    %   07/18 MST - Changed behavior to add a pure "aparc" and "aparc.a2009s" atlas (see Note_2)
    %
    % See also: afni_extract_whereami, atlas_table_fsaverage, atlasFSTranslate, 
    
    if nargin < 3, overwrite = 0; end
    
    MAX_DIST = 10;               % Will report as Not Found if no atlas within 5 mm
    USE_DURAL_PROJ_BACKUP = 0;  % TODO: If a channel Not Found within 5 mm, use the dural projection as backup and check again

    
    element_info = getElementInfo(subj, rootEEGdir);
   
    locDirs = localizer_create_directories(subj, rootEEGdir);
    fname_monopolar = fullfile(locDirs.atlas, 'atlas_monopolar.csv');
    if exist(fname_monopolar, 'file') && ~overwrite
        fprintf('Atlas file found (monopolar), will not overwrite: %s\n', fname_monopolar);
        
    else
        fprintf('Writing atlas file (2+ minutes...)\n');
    
        fname_leads = fullfile(locDirs.mr_pre, 'leads.csv');
        assert(exist(fname_leads,'file') > 0, 'File not found: %s\n', fname_leads);
        leads = readtable(fname_leads);

        xyz     = leads{:,{'x','y','z'}};
        chans   = leads.chanName;
        atlas_table = afni_extract_whereami(...
            'working_dir',  locDirs.working, ...
            'suma_dir',     fullfile(locDirs.fs_subj, 'SUMA'), ...
            'coords',       xyz, ...
            'coord_name',   chans, ...
            'orientation',  'RAI', ...
            'overwrite',    1);
        if isempty(atlas_table)
            error('Atlas error: Problem using AFNI whereami command');
            % Note: Mike has encountered a problem with afni not registering the ORIG template space
            %   Check if this is your problem. This particular problem can be fixed by adding this block
            %   to the SessionAtlases.niml file or your ~/abin/AFNI_atlas_spaces.niml:
            % <TEMPLATE_SPACE
            %   space_name="ORIG"
            %   generic_space="ORIG"
            %   comment="Native space of subject"
            % ></TEMPLATE_SPACE>
        end

        % Come up with the list of all aparc (Desikan) and aparc.a2009s (Destrieux) labels
        [~,desikan] = isInDKAtlas();
        [~,destrieux] = isInDestAtlas();
        parc_only_labls = union(desikan, destrieux);
        
        
        [new_atlas_table, not_found_chan] = process_table(atlas_table, chans, subj, element_info, MAX_DIST, parc_only_labls);
        
        % Rerun extraction with just those channels not found using a snapped-to-dural coords (if possible)
        if USE_DURAL_PROJ_BACKUP && ~isempty(not_found_chan)
            fname = fullfile(locDirs.mr_pre, 'coords_snap_1.csv');
            if exist(fname, 'file')
                snap = readtable(fname);
                snap = snap(ismember(snap.chanName, not_found_chan), :);
                atlas_table = afni_extract_whereami(...
                    'working_dir',  locDirs.working, ...
                    'suma_dir',     fullfile(locDirs.fs_subj, 'SUMA'), ...
                    'coords',       snap{:, {'x','y','z'}}, ...
                    'coord_name',   snap.chanName, ...
                    'orientation',  'RAI', ...
                    'overwrite',    1);
                [refound_atlas_table, not_found_chan] = process_table(atlas_table, chans, subj, element_info, MAX_DIST);
            end
        end

        for i = 1:numel(not_found_chan)
            new_atlas_table = addBlankTableRow(new_atlas_table);
            new_atlas_table(end,:).chanName = not_found_chan(i);
        end

        % save
        writetableSafe(new_atlas_table, fname_monopolar);
        fprintf('written!\n');
    end
    
 
    % do again for bipolar xyz
    fname_biplar = fullfile(locDirs.atlas, 'atlas_bipolar.csv');
    fname_mid_euc = fullfile(locDirs.mr_pre, 'coords_mid_euclid.csv');
    if exist(fname_biplar, 'file') && ~overwrite
        fprintf('Atlas file found (bipolar), will not overwrite: %s\n', fname_biplar);
    else
        if exist(fname_mid_euc, 'file')
            
            fprintf('Writing bipolar atlas file (2+ minutes...)');
            leads = readtable(fname_mid_euc);

            xyz     = leads{:,{'x','y','z'}};
            chans   = leads.chanName;
            atlas_table = afni_extract_whereami(...
                'working_dir',  locDirs.working, ...
                'suma_dir',     fullfile(locDirs.fs_subj, 'SUMA'), ...
                'coords',       xyz, ...
                'coord_name',   chans, ...
                'orientation',  'RAI', ...
                'overwrite',    1);
            
            % Come up with the list of all aparc (Desikan) and aparc.a2009s (Destrieux) labels
            [~,desikan] = isInDKAtlas();
            [~,destrieux] = isInDestAtlas();
            parc_only_labls = union(desikan, destrieux);
        
            [new_atlas_table, not_found_chan] = process_table(atlas_table, chans, subj, element_info, MAX_DIST, parc_only_labls);

            for i = 1:numel(not_found_chan)
                new_atlas_table = addBlankTableRow(new_atlas_table);
                new_atlas_table(end,:).chanName = not_found_chan(i);
            end
            
            % save
            writetableSafe(new_atlas_table, fname_biplar);
            fprintf('written!\n');
        
        else
            fprintf('Cannot atlasize bipolars; file not found: %s\n', fname_mid_euc);
        end
    end
end

function [t,not_found] = process_table(atlas_table, chans, subj, element_info, MAX_DIST, parc_only_labels)
    
    
    
    % the table will have fields: atlas, within_mm, label, (prob_type, code), coord_name, 
    atlases = {'aseg','aparc+aseg','aparc.a2009s+aseg'}; % volume (DK), Desikan-Killiany, and Destrieaux, respectively
    atlases = [atlases, {'aparc', 'aparc.a2009s'}]; % MST 07/18 - add these derived, non-segmentation atlases
    atlas_table = atlas_table(atlas_table.within_mm <= MAX_DIST, {'atlas','within_mm','label', 'coord_name', 'code'});
    
    t       = [];
    rows    = [];
    not_found = [];
    
    % for every channel
    for c = 1:numel(chans)

        tagName = char(util_split_stringnum(chans{c}));
        if contains(tagName, '-')
            temp = strsplit(tagName,'-');
            tagName = temp{1};
        end
        hwType = char(element_info(strcmpi(element_info.tagName, tagName), :).hardwareType);
        
        
        chan_table = atlas_table(strcmpi(atlas_table.coord_name, chans{c}), :);
        if isempty(chan_table)
            warning('Channel not found: %s', chans{c});
            not_found = [not_found chans(c)];
            continue
        end
        
        % find the nearest lobe
        t_temp = sortrows(chan_table, 'within_mm');
        desikan_lobe = cellfun(@aparc2lobe, t_temp.label, 'uniformOutput',0);
        % find first row with this channel that had a desikan-killiany lobe
        r = find(~strcmpi(desikan_lobe, 'unknown'), 1);
        if r > 0
            lobe = desikan_lobe{r};
        else
            lobe = '';
        end
        
        % for each of the atlases
        for i = 1:numel(atlases)
            atlas = atlases{i};
            
            % for our customized parc-only ones, add the expected +aseg for lookup purposes
            if any(cellfun(@(x) strcmpi(atlas, x), {'aparc','aparc.a2009s'}))
                atlasParcOnly = atlas;
                atlas = [atlas, '+aseg'];
                isParcOnly = 1;
            else
                atlasParcOnly = '';
                isParcOnly = 0;
            end
            
            atlas_chan_table = chan_table(strcmpi(chan_table.atlas, atlas),:);       
            
            % some subjects (probably older runs) only have _rank atlases,
            % the only difference is that the code has been renumbered from FS by afni 
            % Right now the code isn't used, so this is (for the moment) okay:
            if isempty(atlas_chan_table)
                atlas = [atlas '_rank'];
                atlas_chan_table = chan_table(strcmpi(chan_table.atlas, atlas),:);
            end
            
            % Filter out unknown labels
            atlas_chan_table = atlas_chan_table(~contains(lower(atlas_chan_table.label), 'unknown'), :);
            
            % Filter out non-parc (ie seg) labels
            if isParcOnly
                % NOTE: labels at this point may have prefixes like "ctx-lh-", but parc_only_labels don't
                %       have these prefixes. Using the "contains" function accounts for these partial matches:
                parcOnly_mask = contains(atlas_chan_table.label, parc_only_labels, 'ignoreCase',1);
                atlas_chan_table = atlas_chan_table(parcOnly_mask, :);
                atlas_chan_table.atlas = repmat(cellstr(atlasParcOnly), height(atlas_chan_table), 1);
            end
            
            % only report the nearest 1 location nearest to this channel for each atlas
            m = min(atlas_chan_table.within_mm);
            atlas_chan_table = atlas_chan_table(atlas_chan_table.within_mm == m, :);
            
            % it is possible that there are >1 equally near locations. For each of them...
            for k = 1:height(atlas_chan_table)
                row = table2struct(atlas_chan_table(k, :));
                
                parts = [];
                s = struct();
                s.chanName  = row.coord_name;
                s.hemisphere= '';
                s.lobe      = lobe;
                s.dist_mm   = row.within_mm;
                s.label     = row.label;
                s.hardwareType = hwType;
                s.class     = '';
                s.label_trim= row.label;
                s.desc      = row.label;
                s.code      = row.code;
                s.atlas     = row.atlas;

                if contains(s.label, 'Right') || contains(s.label, 'rh')
                    s.hemisphere = 'rh';
                    s.label_trim = strrep(s.label_trim, 'Right', '');
                elseif contains(s.label, 'Left') || contains(s.label, 'lh')
                    s.hemisphere = 'lh';
                    s.label_trim = strrep(s.label_trim, 'Left', '');
                end
                
                switch s.atlas
                    case {'aseg' 'aseg_rank'}
                        
                    case {'aparc+aseg' 'aparc.lobesStrict+aseg' 'aparc+aseg_rank' 'aparc'}
                        parts = strsplit(s.label,'-');
                        
                    case {'aparc.a2009s+aseg' 'aparc.a2009s+aseg_rank', 'aparc.a2009s'}
                        parts = strsplit(s.label,'_');
                        s.desc = atlasFSTranslate('destrieux', s.label, 'silent',1);
                end
                
                if ~isempty(parts)
                    if strcmp(parts{1}, 'ctx')
                        parts = parts(2:end);
                    elseif strcmp(parts{1}, 'wm')
                        parts = parts(2:end);
                    end
                    
                    if strcmp(parts{1}, 'lh')
                        s.hemisphere = 'lh';
                        parts = parts(2:end);
                    elseif strcmp(parts{1}, 'rh')
                        s.hemisphere = 'rh';
                        parts = parts(2:end);
                    end 
                    if strcmp(parts{1}, 'G') && numel(parts) > 2 && strcmp(parts{3}, 'S')
                        s.class = 'gyrus-sulcus';
                        parts = parts(4:end);
                    elseif strcmp(parts{1}, 'S')
                        s.class = 'sulcus';
                        parts = parts(2:end);
                    elseif strcmp(parts{1}, 'G')
                        s.class = 'gyrus';
                        parts = parts(2:end);
                    end
                    s.label_trim = strrep(strjoin(parts,'-'), '_', '-');
                end
                s.label_trim = strip(strip(s.label_trim, 'both','-'), 'both','-');
                
                switch s.atlas
                    case {'aseg' 'aseg_rank'}
                        
                    case {'aparc+aseg' 'aparc.lobesStrict+aseg' 'aparc+aseg_rank' 'aparc'}
                        s.desc = atlasFSTranslate('desikan', s.label_trim, 'silent',1);
                        
                    case {'aparc.a2009s+aseg' 'aparc.a2009s+aseg_rank', 'aparc.a2009s'}
                        
                end
                
                
                if isempty(rows)
                    rows = s;
                else
                    rows(length(rows)+1) = s;
                end
            end % end of atlas_chan per-row loop
        end % end of atlas
        
        
    end % channel loop
    
    t = struct2table(rows);

    
    % 1) not sure why there are row duplicates
    % 2) unique throws an error on the desc and hardware type fields sometimes. Therefore:
    vars = setdiff(t.Properties.VariableNames, {'desc', 'hardwareType'});
    [~, it] = unique(t(:,vars), 'stable');
    t = t(it, :);
    
end

% If we want to ever add the afni-independent freesurfer-based surface parcellation, base it off of this code:
% bd = braindata('suma_dir', fullfile(locDirs.fs_subj, 'SUMA'));
%         parc_table = [];
%         seg_table = table();
%         for hemisphere = hems
%             hemi_chans_subd     = getLeads(subj, root, 'whichHemi',[hemisphere 'h'], 'hardwareType',{'subdural','micro-subdural'}, 'jacktable',jacktable);
%             hemi_chans_depth    = getLeads(subj, root, 'whichHemi',[hemisphere 'h'], 'hardwareType',{'depth'}, 'jacktable',jacktable);
%             bd.fn_load_surf('hemi',[hemisphere 'h'], 'surfType','pial');
%             bd.fn_load_parc('parcType','a2009s','hemi',[hemisphere 'h']);
%             lut                 = bd.parc.a2009s.lut;
%             lead_roi_lut_hemi   = lead_roi_lut(ismember(lead_roi_lut.chanName, hemi_chans_subd),:);
%             parc_table_hem      = lead_roi_lut_hemi(:, {'chanName'});
% 
%             % Cocjin magic;
%             ind = table();
%             [ind.label ind.node] = util_group2ind(lut.label_name, lut.vert_idx); % break out into individual format
%             parc_table_hem.parc = cellfun(@(nodes) unique(ind.label(ismember(ind.node, nodes))), lead_roi_lut_hemi.lead_verts, 'UniformOutput',false);
%             parc_table = [parc_table; parc_table_hem];
% 
% 
%         end