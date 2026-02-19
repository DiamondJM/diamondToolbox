function toolbox_dir = getToolboxDir()
    % getToolboxDir returns the location of the user's eeg_toolbox folder
    
    thispath = mfilename('fullpath');
    assert(contains(thispath, 'eeg_toolbox'), 'This function should be a descendent of the eeg_toolbox directory');
    
    isTrunk = contains(thispath, 'trunk');
    
    % get parent dir
    
    isTrunk = contains(thispath, 'trunk');
    
    while contains(thispath, 'eeg_toolbox')
        thispath = fileparts(thispath);
    end
    
    if isTrunk
        toolbox_dir = fullfile(thispath, 'eeg_toolbox', 'trunk');
    else
        toolbox_dir = fullfile(thispath, 'eeg_toolbox');
    end
    
    
end