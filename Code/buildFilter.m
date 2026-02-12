function myFilter = buildFilter(Fs,filterWindow, filterOrder)

warning('Should you be using this? Should only be used for programmatic filter design if one of the pre-saved filters is insufficient.'); 
% We should phase this out, and begin using saved filters designed by
% filterDesigner. 
% The filter should be designed and saved, not built on the fly. 

if nargin == 2; filterOrder = 1999; end
% if nargin == 2; filterOrder = 49; end
% warning('Changed filter order to 49'); 

nyquist = Fs / 2;
transitionWidth = .1;
% transitionWidth = .2; warning('Using transition width %f',transitionWidth);
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
    clf
    
    freqVals = abs(fft(myFilter));
    freqDomain = linspace(0, nyquist,floor((filterOrder + 1) / 2) + 1);
    plot(freqDomain,freqVals(1:length(freqDomain)))
    
    set(gca,'xlim',[0 (1 + transitionWidth)*(filterWindow(2))*1.4])
    
end




end