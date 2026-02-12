function [rh,lh,failureFlag] = retrieveResectedVertex(subj)

[fn, failureFlag] = locateVertices(subj);

if ~failureFlag 
    % rh = readtable(fullfile(fn,'RH_node_to_resection.txt')); rh = logical(rh.Var2); 
    % lh = readtable(fullfile(fn,'LH_node_to_resection.txt')); lh = logical(lh.Var2); 
    
    rh = load(fullfile(fn,'std.141.rh.smoothed.resection.msk.1D.dset')); rh = logical(rh(:,2)); 
    lh = load(fullfile(fn,'std.141.lh.smoothed.resection.msk.1D.dset')); lh = logical(lh(:,2)); 

    if ~(any(lh) || any(rh)); failureFlag = true; end
else
    warning('Resected vertices for %s not found.',subj); 
    lh = []; rh = [];
    return;
end

end