
function dataMaster = extractPhasePow(sourceLoc)

TimeSeries = sourceLoc.TimeSeries;

% Use the following lines to shuffle the lead labels. 
% warning('Shuffling lead labels.') 
% TimeSeries = TimeSeries(:,randperm(size(TimeSeries,2)));

Fs = sourceLoc.Fs;
filterWindow = sourceLoc.paramStruct.filterWindow;

[a,b] = size(TimeSeries); 

if a < b
    warning('Time series seems to be wrong dimensions. Transposing. \n');
    TimeSeries = TimeSeries';
end


myFilter = buildFilter(Fs,filterWindow,1999);
filtData = filtfilt(myFilter,1,TimeSeries);
hilbertData = hilbert(filtData);
phase = angle(hilbertData);

freqData = Fs * diff(unwrap(phase)) / (2 * pi);
freqData = medfilt1(freqData,1000,'truncate');
freqData(end + 1,:) = freqData(end,:);

%% Process and pack 

pow = abs(hilbertData);
realData = real(filtData);

dataMaster.phase = phase;
dataMaster.pow = pow;
dataMaster.freq = freqData;
dataMaster.real = realData; 



