function myCI = calculateConfInterval(y)

y(isnan(y)) = []; 

if min(size(y)) == 1; N = length(y); 
else; N = size(y,1); 
end
% yMean = mean(y,'omitnan');
ySem = std(y,'omitnan') / sqrt(N); 

ts = tinv(0.975,N - 1); 
myCI = ts * ySem;

end
