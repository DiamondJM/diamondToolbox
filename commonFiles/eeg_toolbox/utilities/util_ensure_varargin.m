function my_var = util_ensure_varargin(my_var)
% Desc: ensure that your varargin is in proper name-value pair form
%       if you want to input it as an info struct, you can!
%       (which is actually really nice)
% Inputs:
%   my_var: cell array of either name-value pair or of a single info struct
% Outputs:
%   my_var: varargin as a cell array of name-value pair form
% Notes:
%   it's kinda silly...
%   in practice... this will convert a struct to a cell array BACK to a
%   struct
%   really so default inputs can be used properly...
%%
if ~isempty(my_var)
    switch class(my_var{1})
        case 'string'
        case 'struct'
            my_cnt = 0;
            new_var = {};
            for my_field = fields(my_var{1})'
                my_cnt = my_cnt + 1;
                new_var{1,my_cnt} = my_field{:};
                
                my_cnt = my_cnt + 1;
                new_var{1,my_cnt} = my_var{1}.(my_field{:});
            end
        case 'cell'
            if isempty(my_var{1})
                my_var = {};
            end
    end
end
%%
return