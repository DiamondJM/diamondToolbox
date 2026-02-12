function inputs = util_extract_inputs(default_cell, input_pair_list)
% Desc: extract inputs into a structure given a default table and
%       input_pair_list
% Inputs:
%   default_table: cell array, contains columns that correspond to...
%                  param: parameter name
%                  default: default value for that parameter
%                  vartype: variable type, 'optional', 'required', or
%                  'parameter'
%   input_pair_list: cell array, name-value style
%                    e.g. {'input', [1], 'output', [2]}
% Outputs:
%   inputs: structure where fields are of 'param' name and either 'default'
%   or 'input_pair_list' value
% Example:
%   fo lata

% make it so we can take in either a struct or an input pair list
input_pair_list = util_ensure_varargin(input_pair_list);

% split out table into defaults
default_table = cell2table(default_cell,'VariableNames',{'param','default','vartype'});

% create default inputParser object
ip = inputParser; ip.KeepUnmatched = 1; % keep unmatched to handle when a parameter isn't in defaults
                                        % may consider stroing the
                                        % unmatched within 'inputs' 
                                        % so debugging errors in typing
                                        % is a lil easier
for i_row = 1:height(default_table)
    my_row = default_table(i_row,:);
    switch my_row.vartype{:}
        case 'optional'
            which_func = 'addOptional';
        case 'parameter'
            which_func = 'addParameter';
        case 'required'
            which_func = 'addRequired';
    end
    if isstruct(my_row.default)
        ip.(which_func)(my_row.param{:}, my_row.default);
    else
        try
            ip.(which_func)(my_row.param{:}, my_row.default{:});
        catch
            ip.(which_func)(my_row.param{:}, my_row.default);
        end
    end
end

% replace stuff that you provided
ip.parse(input_pair_list{:});

% extract these into your inputs structure
for i_row = 1:height(default_table)
    which_param = default_table{i_row,'param'}{:};
    inputs.(which_param) = ip.Results.(which_param);
end

for which_param = fieldnames(ip.Unmatched)'
    inputs.(which_param{:}) = ip.Unmatched.(which_param{:});
end
% inputs = util_flood_list(inputs,inputs);
return