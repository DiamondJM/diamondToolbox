function meanSdTTString = printMeanSdTT(groupA,groupB,varargin)


p = inputParser;
addParameter(p,'usePaired',false)
addParameter(p,'tail','both'); 
parse(p,varargin{:})
usePaired = p.Results.usePaired;
tail = p.Results.tail; 


meanSdA = printMeanSD(groupA);
if exist("groupB",'var') && ~isempty(groupB)
    meanSdB = printMeanSD(groupB); 
    TTstring = printTTest(groupA,groupB,'usePaired',usePaired,'tail',tail);


else
    TTstring = printTTest(groupA,[],'usePaired',usePaired,'tail',tail);



end

%% Print

if ~exist("groupB",'var')
    meanSdTTString = sprintf('$%s$; %s',meanSdA, TTstring); 
else
    meanSdTTString = sprintf('$%s \\text{ versus } %s$; %s',meanSdA,meanSdB, TTstring);
end

disp(meanSdTTString); 