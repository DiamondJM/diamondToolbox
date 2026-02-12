function currentVec = revertToVector(currentVec)


[a,b] = size(currentVec);

if a ~= 1 && b ~= 1; error('Input must be a vector'); end

if a > b; currentVec = currentVec'; end

end