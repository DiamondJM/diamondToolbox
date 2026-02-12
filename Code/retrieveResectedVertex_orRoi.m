function [rh,lh,failureFlag] = retrieveResectedVertex_orRoi(subj)


% Hopefully this can be phased out, and replaced with
% retrieveResectedVertex, once Price's masks are finished. 

[rh,lh,failureFlag] = retrieveResectedVertex(subj);
if ~failureFlag; return; end 

%% Now we can try to do it the old way? 

[~, failureFlag] = retrieveMaskValues(subj);
if failureFlag; return; end 

[~, roi, ~, isLeftResection]  = retrieveResectionRois(subj); 

[~, myBp] = retrieveBraindata(subj); 

[rh,lh] = deal(false(myBp.stdNumVert,1)); 
if isLeftResection; lh(roi) = true; 
else; rh(roi) = true; 
end


end