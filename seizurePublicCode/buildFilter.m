function myFilter = buildFilter(Fs,filterWindow, filterOrder)


if nargin == 2; filterOrder = 1999; end

nyquist = Fs / 2;
transitionWidth = .1;
% transitionWidth = 1;
% warning('Using transition width 1');
% filterWindow = [5 10];

idealResponse = [0 0 1 1 0 0];
fFrequencies = [0 ...
    (1 - transitionWidth) * filterWindow(1) ...
    filterWindow(1)...
    filterWindow(2) ...
    (1 + transitionWidth) * filterWindow(2) ...
    nyquist] / nyquist;
myFilter = firls(filterOrder, fFrequencies,idealResponse);


checkFilter = false;

if checkFilter
    figure
    
    freqVals = abs(fft(myFilter));
    freqDomain = linspace(0, nyquist,floor((filterOrder + 1) / 2) + 1);
    plot(freqDomain,freqVals(1:length(freqDomain)))
    
    set(gca,'xlim',[0 (1 + transitionWidth)*(filterWindow(2))*1.4])
    
end




end