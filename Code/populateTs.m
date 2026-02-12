function mySingleton = populateTs(mySingleton,varargin)

% General populateTs for depths and subdurals. 

% Comments
% Primarily, the point of this code is to divide the given structure into
% depths and subdurals. 
% If leads exist, Ts is attempted to be loaded. If the Ts doesn't exist,
% it's created afresh. Importantly, regardless of whether the existing Ts
% is loaded or created afresh, the function will TRY to match the given
% leads. 
% If the leads don't match, it's created anew. 
% Sometimes the leads CAN'T match. In this case, a message will be given,
% and this just means the struct should be re-processed, ideally, just pull
% event from server and run eventManager on it. 

% If the leads don't exist, the Ts and leads are loaded from whatever is
% saved in the directory. 

% If forceTs is passed, as of now, this will wipe prior channels. I think
% this is the best behavior. 
% To make it otherwise, and try to preserve channel ordering and just rerun
% the time series, refer to populateTs_byHardwareType. p = inputParser;

%% Preamble 
assert(length(mySingleton)==1,'Should be a singleton.'); 
assert(isequal(mySingleton.task,'seizure')); 

p = inputParser;
addParameter(p,'useBaseline',false);
addParameter(p,'forceTs',false);
parse(p,varargin{:})
forceTs = p.Results.forceTs;
useBaseline = p.Results.useBaseline;

[~,hasSubdural,hasDepth,leadsSubdural,leadsDepth] = generateChansMasterList(mySingleton.subjID);

if ~isfield(mySingleton,'sourceLoc'); mySingleton.sourceLoc = struct; end
if isfield(mySingleton.sourceLoc,'TimeSeries') && ~forceTs && ~useBaseline; return; end
if isfield(mySingleton.sourceLoc,'chanNames'); chansPassed = true; chanNamesPassed = mySingleton.sourceLoc.chanNames; 
else; chansPassed = false; 
end

%% Split

if hasSubdural
    mySingletonSubdural = mySingleton;
    if chansPassed
        subduralInds = ismember(chanNamesPassed,leadsSubdural);
        mySingletonSubdural.sourceLoc.chanNames = mySingletonSubdural.sourceLoc.chanNames(subduralInds);
    end
    mySingletonSubdural = populateTs_byHardwareType(mySingletonSubdural,'hardwareType','subdural','forceTs',forceTs,'useBaseline',useBaseline);
end
if hasDepth
    mySingletonDepth = mySingleton;
    if chansPassed
        depthInds = ismember(chanNamesPassed,leadsDepth);
        mySingletonDepth.sourceLoc.chanNames = mySingletonDepth.sourceLoc.chanNames(depthInds);
    end
    mySingletonDepth = populateTs_byHardwareType(mySingletonDepth,'hardwareType','depth','forceTs',forceTs,'useBaseline',useBaseline);
end

%% Combine 

if hasSubdural && hasDepth
    % If both exist, we keep the subdural.
    mySingletonSubdural.sourceLoc.TimeSeries = [mySingletonSubdural.sourceLoc.TimeSeries mySingletonDepth.sourceLoc.TimeSeries];
    mySingletonSubdural.sourceLoc.chanNames = [mySingletonSubdural.sourceLoc.chanNames; mySingletonDepth.sourceLoc.chanNames];
    assert(isequal(mySingletonSubdural.sourceLoc.Fs,mySingletonDepth.sourceLoc.Fs));
    mySingleton = mySingletonSubdural;
elseif hasSubdural
    mySingleton = mySingletonSubdural;
elseif hasDepth
    mySingleton = mySingletonDepth;
end
    
end