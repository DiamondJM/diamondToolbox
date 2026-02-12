
function [roic, roi, resectedVertex, isLeftResection]  = retrieveResectionRois(subj) 

warning('off','stats:pdist2:DataConversion');


warning('This should only be used if Price''s mask not found');
% Should be phased out if at all possible. 
[myBd, myBp] = retrieveBraindata(subj); 

isResected = myBd.docs.jacksheet.isResected;
if isequal(class(myBd.docs.jacksheet.isResected),'cell')
    isResected = str2double(isResected);
end
isResected = logical(isResected);
isResected = myBd.docs.jacksheet.chanName(isResected);
[lResec, rResec] = myBd.chans2hem(isResected, myBd.docs.jacksheet);

if length(lResec) > length(rResec)
    isLeftResection = true;
    resecSurfVert = myBp.surfaces.pial_lh.vertices;
    surfString = 'lh';
else
    isLeftResection = false;
    resecSurfVert = myBp.surfaces.pial_rh.vertices;
    surfString = 'rh';
end

% if exist(sprintf('resectedVertex/%s.mat',subj),'file')
fName = sprintf('/Users/diamondjm/Documents/Josh_work/commonFiles/resectedVertex/%s.mat',subj);
if exist(fName,'file')
    load(fName,'resectedVertex');
else
    [maskValues, failureFlag] = retrieveMaskValues(subj);
    
    if failureFlag; error('Failed to find resection mask.'); end
    
    k = boundary(maskValues);
    k = unique(k);
    maskValues = maskValues(k,:);
    
    
    
    resectedVertex = pdist2(resecSurfVert,maskValues);
    resectedVertex = find(any(resectedVertex < 0.5,2));
    
    save(fName,'resectedVertex');
end
    
[roic,roi]  = myBd.vertex2ROI(resectedVertex, surfString);


end
