function [d1, d0] = coeffToOrientation(coeff)


% Orientation is RAS
oPos = {'Right','Anterior','Superior'};
oNeg = {'Left','Posterior','Inferior'};
d1 = cell(3,3);
d0 = d1;

if size(coeff,2) < 3 
    coeffTemp = zeros(3); 
    coeffTemp(:,1:size(coeff,2)) = coeff;
    coeff = coeffTemp;
end

orientationPC = eye(3) * coeff';

for ii = 1:3
    [~,inds] = sort(abs(orientationPC(ii,:)),'descend');
    
    posInds = sign(orientationPC(ii,inds)) == 1;

    
    d1(ii, posInds) = oPos(inds(posInds));
    d1(ii, ~posInds) = oNeg(inds(~posInds));
    
    d0(ii, posInds) = oNeg(inds(posInds));
    d0(ii, ~posInds) = oPos(inds(~posInds));
    
    
    
end


% Each ROW of d provides the biggest, medium, and smallest influence of the
% initial dimensions, into that row. 



end