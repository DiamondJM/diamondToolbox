function [f,sMean] = checkLineNoise(x,Fs,plotting) 


if nargin == 2; plotting = true; end 

% if plotting; clf; end
x(isnan(x)) = 0;

[s,f,t] = spectrogram(x,5000,0,[],Fs,'power');
% ind = find(f>128,1,'first');
% ind = find(f>450,1,'first');
ind = length(f); 
f = f(1:ind);
s = s(1:ind,:);
sMean = mean(abs(s),2);

if ~plotting; return; end 

plot(f,sMean)
hold on
plot([60 60],ylim,'k:'); plot([120 120],ylim,'k:'); 

end