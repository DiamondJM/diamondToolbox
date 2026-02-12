function patients = determinePatientsAvailability(varargin)


%% Preamble 

p = inputParser;
% addParameter(p,'requireMask',false);
addParameter(p,'useOutcome',false);
addParameter(p,'useUnion',false); 
parse(p,varargin{:})
% requireMask = p.Results.requireMask;
useOutcome = p.Results.useOutcome;
useUnion = p.Results.useUnion; 


%% Onto the function 

patientsSeizure = determinePatientsFromDirectory(pullSeizureDirectory);
% warning('Ignoring seizure.'); 
patientsSpikes = determinePatientsFromDirectory(pullSpikesDirectory);

if useUnion; patients = union(patientsSeizure,patientsSpikes);
else; patients = intersect(patientsSeizure,patientsSpikes);
end

% patients = patientsSpikes;

validInds = true(size(patients)); 

if useOutcome
    hasOutcome = false(size(patients));
    hasMask = false(size(patients));
    
    for ii = 1:length(patients)
        subj = sprintf('NIH0%d',patients(ii)); 
        [~,failureFlag] = retrieveOutcome(subj);
        hasOutcome(ii) = ~failureFlag;
        
        
        [~,failureFlag] = locateVertices(subj);
        if failureFlag
            [~,failureFlag] = retrieveMaskValues(subj);
        end
        hasMask(ii) = ~failureFlag;
        
    end
    validInds = validInds & hasOutcome; 
    validInds = validInds & hasMask;
end

patients = patients(validInds);


end
