function [ts,mu,sigma] = zScore_omitNan(ts,dim)

if nargin == 1; dim = 1; end 

[a,b] = size(ts); 
if a < b && dim == 1
    if a > 1; error('Check input dimensions.');
    elseif a == 1; dim = 2; 
    end
end


mu = mean(ts,dim,'omitnan'); 
sigma = std(ts,[],dim,'omitnan');

ts = (ts - mu) ./ sigma; 

end
