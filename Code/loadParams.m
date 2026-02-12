function paramStruct = loadParams(params,paramStruct)

% Params is the new set of parameters, loaded onto paramStruct. 

newParams = fieldnames(params);

% oldParams = fieldnames(paramStruct);

% newParams = setdiff(newParams,oldParams);
% The fuck is the point of the above? 
% We should be squashing old params, that is, overwriting them, not letting
% them be. Failure to overwrite was causing problems. 
% But since I haven't checked the above line, let's leave a warning. 
% warning('Overwriting old parameters.'); 

for ii = 1:length(newParams)
    
    paramStruct.(newParams{ii}) = params.(newParams{ii});
    
end
