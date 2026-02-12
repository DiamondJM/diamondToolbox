function pcaInfo = generatePcaInfo(locationsMaster,varargin)

%% 

p = inputParser;
addParameter(p,'useEye',false);
parse(p,varargin{:})
useEye = p.Results.useEye;


%%

if size(unique(locationsMaster,'rows'),1) <= 2     
    warning('Insufficient points to perform PCA.');
    useEye = true; 
end

if useEye
    coeff = eye(3);
    mu = zeros(1,3);
else
    [coeff, ~, ~, ~, ~, mu] = pca(locationsMaster,'centered',true);
end

% coeff(:,[1 2 3]) = coeff(:,[2 1 3]);
% To flip x and y 

pcaInfo = struct;
pcaInfo.coeff = coeff;
pcaInfo.mu = mu; 

end