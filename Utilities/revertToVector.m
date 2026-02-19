function currentVec = revertToVector(currentVec)


[a,b] = size(currentVec);

if a ~= 1 && b ~= 1
    warning('Input must be a vector');currentVec = nan; 
end

if a > b; currentVec = currentVec'; end

end