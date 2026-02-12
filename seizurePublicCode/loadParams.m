function paramStruct = loadParams(params,paramStruct)

newParams = fieldnames(params);

oldParams = fieldnames(paramStruct);

newParams = setdiff(newParams,oldParams);

for ii = 1:length(newParams)
    
    paramStruct.(newParams{ii}) = params.(newParams{ii});
    
end
