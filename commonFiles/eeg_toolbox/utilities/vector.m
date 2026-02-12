function V = vector(A)
    % returns A reshaped to a column vector
    sz = numel(A);
    V = reshape(A,sz,1);
end