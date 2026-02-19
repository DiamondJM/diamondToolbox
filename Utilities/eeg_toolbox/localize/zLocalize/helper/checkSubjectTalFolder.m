function [success, missing_anything, needs_deleting] = checkSubjectTalFolder(subj, root, do_delete)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin < 1
        subj = 'NIH044';
    end
    if nargin < 2 
        dirs = evalin('base','dirs');
        root = dirs.zloc;
    end
    if nargin < 3
        do_delete = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Z


    tal  = fullfile(root, subj, 'tal');
    zloc = fullfile(tal, 'zloc');

    zloc_dirs_need = { ...
                    'anchors'
                    'CT_1'
                    'CT_2'
                    'freesurfer'
                    'logs'
                    'metrics'
                    'mr_pre'
                    'working'};

    zloc_dirs_ok_prefix = { ...
                      'ALICE'
                      'mr'};

    tal_files_need = { ...
                      'atlas/atlas_bipolar.csv'
                      'atlas/atlas_monopolar.csv'
                      'roi/lead_ROIC_LUT_monopolar.mat'
                      'roi/lead_ROIC_LUT_bipolar.mat'
                      'roi/lead_mesh_LUT_monopolar.mat'
                      'roi/lead_mesh_LUT_bipolar.mat'
                      'roi/ROIC_mesh_LUT_*.csv'};

    tal_files_bad = { ...
                    'roi/ROIC_lead_LUT.mat'
                    'roi/Mesh_LUT.csv'
                    'roi/ROIC_Mesh_LUT.csv'
                    'roi/ROIC_lead_LUT_bipolar.mat'
                    'roi/ROIC_lead_LUT_monopolar.mat'};
                  
                  
    zloc_files_bad = { ...
                    'anchors/coords.mr_pre.anchors.csv'
                    'anchors/*.gii'
                    'anchors/*.nii'
                    'roi*'
                    'atlas*'
                    'working/temp_*'
                    'CT_1/img_check*'
                    'CT_2/img_check*'};

    lsCount = @(f) numel(lsCell(f));
    lsAny   = @(f) lsCount(f) > 0;
    lsAnyCell = @(c) cellfun(lsAny, c);
    isdirCell = @(c) cellfun(@isdir, c);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Z
    fprintf('\n===================================== %s\n', subj);

    needs_deleting = 0;
    missing_anything = 0;



    % Dirs necessary
    found = isdirCell(fullfile(zloc, zloc_dirs_need));
    if ~all(found)
        fprintf('The following %s zloc directories were missing:\n', subj);
        disp(zloc_dirs_need(~found));
        missing_anything = 1;
    end

    % Dirs that need deleting
    present = lsCell(zloc);
    extra = setdiff(present, zloc_dirs_need);
    for i=1:numel(zloc_dirs_ok_prefix)
        extra = extra(~strncmp(extra, zloc_dirs_ok_prefix{i}, length(zloc_dirs_ok_prefix{i})));
    end
    if ~isempty(extra)
        fprintf('The following %s zloc directories need to be deleted:\n', subj);
        disp(extra(:));
        if do_delete
            for i=1:numel(extra)
                try
                    rmdir(fullfile(zloc, extra{i}));
                    fprintf('\t%s deleted\n', extra{i});
                catch
                    if ~isempty(lsCell(fullfile(zloc, extra{i})))
                        fprintf('Directory %s has contents: \n', fullfile(zloc, extra{i}));
                        disp(lsCell(fullfile(zloc, extra{i})));
                        
                        DELETE_EVERYTHING = 0; % <-- IF 1, NO PROMPT. If 0, ask/warn user
                        if DELETE_EVERYTHING || inputYN('Delete everything?')
                            rmdir(fullfile(zloc, extra{i}),'s');
                            fprintf('\t%s deleted\n', extra{i});
                        else
                            needs_deleting = 1;
                        end
                    end
                end

            end
        else
            needs_deleting = 1;
        end
    end

    %% Files that need deleting
    found = find(lsAnyCell(fullfile(zloc, zloc_files_bad)));
    if ~isempty(found)
        fprintf('The following %s zloc/ files need to be deleted:\n', subj);
        for i=1:numel(found)
            name = fullfile(zloc, zloc_files_bad{found(i)});
            files = lsCell(name);
            disp(files(:));
            if do_delete
                for j=1:numel(files)
                    if isdir(fullfile(fileparts(name), files{j}))
                        rmdir(fullfile(fileparts(name), files{j}), 's');
                    else
                        delete(fullfile(fileparts(name), files{j}));
                    end
                    fprintf('\t%s deleted\n', files{j});
                end
            else
                needs_deleting = 1;
            end
        end
    end
    found = find(lsAnyCell(fullfile(tal, tal_files_bad)));
    if ~isempty(found)
        fprintf('The following %s tal/ files need to be deleted:\n', subj);
        for i=1:numel(found)
            name = fullfile(tal, tal_files_bad{found(i)});
            files = lsCell(name);
            disp(files(:));
            if do_delete
                for j=1:numel(files)
                     if isdir(fullfile(fileparts(name), files{j}))
                        rmdir(fullfile(fileparts(name), files{j}), 's');
                    else
                        delete(fullfile(fileparts(name), files{j}));
                    end
                    fprintf('\t%s deleted\n', files{j});
                end
            else
                needs_deleting = 1;
            end
        end
    end

    % Files necessary
    found = lsAnyCell(fullfile(tal, tal_files_need));
    if ~all(found)
        fprintf('The following %s tal files were missing:\n', subj);
        disp(tal_files_need(~found));
        missing_anything = 1;
    end

    
    success = ~missing_anything && ~needs_deleting;

end