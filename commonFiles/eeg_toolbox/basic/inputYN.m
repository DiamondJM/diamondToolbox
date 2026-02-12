function yes = inputYN(prompt)
    % prompts for a yes/no
    while true
        response = input([prompt ' [y/n]: '], 's');
        if strncmpi(response, 'yes', min(length(response), length('yes')))
            yes = true;
            return;
        end
        if strncmpi(response, 'no', min(length(response), length('no')))
            yes = false;
            return;
        end
    end
end % funciton inputYN