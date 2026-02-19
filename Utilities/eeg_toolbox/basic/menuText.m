function [nums_selected, strs_selected] = menuText(title, varargin)
% MENUTEXT provides a multi-choice text menu to the user to select from
%
% [nums_selected, strs_selected] = menuText(title, option_1, ..., option_n)
%
% Inputs:
%   title - header for menu
%   items - choices
%
% Optional key-value Input:
%   multiSelect - True (default) or False
%
% toggle choice implementation

    multiSelect = 1;

    if nargin < 2, error('No choice items'); 
    elseif nargin >= 2 && iscell(varargin{1})
        % choices were passed in a cell array
        [nums_selected, strs_selected] = menuText(title, varargin{1}{:}, varargin{2:end});
        return;
    end
    if any(strcmpi(varargin, 'multiSelect'))
        val_ndx = find(strcmpi(varargin, 'multiSelect')) + 1;
        try 
            multiSelect = varargin{val_ndx} + 0;
            varargin = varargin([1:val_ndx-2, val_ndx+1:end]); % valid 2nd param
        catch 
            warning('No interpretable value given for key-item pair multiSelect; assuming you did not pass a value');
            varargin = varargin([1:val_ndx-2, val_ndx:end]);
        end 
            
        
    end
    
    % variable setup
    items = varargin;
    strs_selected = [];
    nums_selected = [];
    userInput = '';
    
    if numel(items) == 1
        nums_selected = 1;
        strs_selected = items{1};
        return;
    end
    
    % loop until (q)uit given
    while isempty(userInput) || ~isnumeric(userInput) || multiSelect
        disp(['+' repmat('-',1,10+length(title)) '+']);
        disp(['|     ',title,'     |'])
        disp(['+' repmat('-',1,10+length(title)) '+']);
        for i = 1 : length(items)
            isSel = ismember(i, nums_selected); 
            mark = turnary(isSel, '*', ' ');
            fprintf(' %c %d - %s\n', mark, i, items{i});
        end
        if multiSelect
            fprintf('   Q - Quit\n\n');
            fprintf('Number of selections made: %d\n\n', length(nums_selected));
            fprintf('\nEnter number or [array] of numbers to toggle choices.\n');
            fprintf('Finished? Enter (Q)uit to continue/exit.\n\n');
        else 
            fprintf('\nPlease enter a *number* to make your selection.\n\n');
        end
        
        
        
        
        userInput = input('Selection: ', 's');
        try userInput = eval(userInput); 
        catch
            if ~strncmpi(userInput, 'QUIT', length(userInput))
                fprintf('Unrecognized matlab expression. Use a real matlab expression!\n')
            end
        end
        
        if isnumeric(userInput) && all(arrayfun(@(x) x <= length(items), userInput))
            if multiSelect
                % treat input as a toggle
                uni = union(nums_selected, userInput);
                int = intersect(nums_selected, userInput);
                nums_selected = setdiff(uni, int);
                nums_selected = intersect(nums_selected, 1:length(items));
                strs_selected = items(nums_selected);
            else
                nums_selected = userInput;
                strs_selected = items{nums_selected};
            end
        end
        
        if multiSelect && strncmpi(userInput, 'QUIT', length(userInput))
            break;
        end
    end
    
    
    

end