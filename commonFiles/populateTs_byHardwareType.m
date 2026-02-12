function mySingleton = populateTs_byHardwareType(mySingleton,varargin)

%% Preamble

p = inputParser;
addParameter(p,'forceTs',false);
addParameter(p,'hardwareType','subdural');
addParameter(p,'useBaseline',false);
parse(p,varargin{:})
forceTs = p.Results.forceTs;
hardwareType = p.Results.hardwareType;
useBaseline = p.Results.useBaseline;
assert(isfield(mySingleton,'SeizureType'),'Must be a seizure.');

if ~isfield(mySingleton,'sourceLoc'); mySingleton.sourceLoc = struct; end
if isfield(mySingleton.sourceLoc,'TimeSeries') && ~forceTs && ~useBaseline; return; end


identifier = sessionToIdentifier(mySingleton);
identifier = sprintf('%s.mat',identifier);

if isfield(mySingleton.sourceLoc,'chanNames') ...
        && isfield(mySingleton.sourceLoc,'Fs') ...
        && ~isempty(mySingleton.sourceLoc.chanNames) ...
        && ~forceTs
    % Optional:
    % The above will prohibit the attempt to match chanNames, if
    % forceTs is passed.
    % In fact I'll leave it on at this point.
    
    chansPassed = true; % chanNames and Fs passed
    chanNamesPassed = mySingleton.sourceLoc.chanNames;
    FsPassed = mySingleton.sourceLoc.Fs;
    
else; chansPassed = false;
end

%% Find old?

if useBaseline; filename = '/Users/diamondjm/Documents/Josh_work/commonFiles/Anne/tsBaseline';
else; filename = '/Users/diamondjm/Documents/Josh_work/commonFiles/Anne/tsSeizure';
end
if isequal(hardwareType,'depth'); filename = [filename 'Depth']; end
if ~isfolder(filename); mkdir(filename); end 
filename = fullfile(filename,identifier);

if exist(filename,'file') == 2 && ~forceTs
    load(filename,'TimeSeries','chanNames','Fs')
    
    if chansPassed
        if ~isequal(chanNames,chanNamesPassed); forceTs = true; end
        if ~isequal(Fs,FsPassed); forceTs = true; end
    end
    
    if ~forceTs
        mySingleton.sourceLoc.TimeSeries = TimeSeries;
        mySingleton.sourceLoc.chanNames = chanNames;
        mySingleton.sourceLoc.Fs = Fs;
    end
    
else
    forceTs = true;
end

%% Create new

if forceTs
    tempPwd = pwd; 
    cd(pullSeizureDirectory); 
    mySingleton = completeTsSeizure(mySingleton,'hardwareType',hardwareType,'useBaseline',useBaseline);
    
    TimeSeries = mySingleton.sourceLoc.TimeSeries;
    chanNames = mySingleton.sourceLoc.chanNames;
    Fs = mySingleton.sourceLoc.Fs;
    
    if isempty(mySingleton.sourceLoc.sess_info)
        % This happens when the file is bad or corrupt.
        % return 
        % Do nothing 
    end
    
    if chansPassed
        % We need to rectify these files then.
        [~,inds] = ismember(chanNamesPassed,chanNames);
        
        if any(~logical(inds))
            missingChannels = setdiff(chanNamesPassed,chanNames);
            fprintf('Channels are missing, including:');
            for jj = 1:length(missingChannels); fprintf(' %s',missingChannels{jj}); end
            fprintf('. mySingleton chanNames has been modified. \n');
            inds = inds(logical(inds));
        else
            assert(isequal(chanNames(inds),chanNamesPassed));
        end
        
        chanNames = chanNames(inds);
        TimeSeries = TimeSeries(:,inds);
        
        % Update the singleton
        mySingleton.sourceLoc.TimeSeries = TimeSeries;
        mySingleton.sourceLoc.chanNames = chanNames;
        mySingleton.sourceLoc.Fs = Fs;
        
        % It's unusual for us to get here.  
    end
    save(filename,'TimeSeries','chanNames','Fs')
    
    cd(tempPwd); 
    
end

end

