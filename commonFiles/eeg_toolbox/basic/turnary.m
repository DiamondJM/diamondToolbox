function output = turnary(condition, resultA, resultB)
    % TURNARY evaluates a condition and returns resultA if true, otherwise result B

    % Return true conditions
    iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
    
    if condition == true
        output = resultA; 
    else
        output = resultB;
    end

end