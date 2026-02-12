function failureFlag = saveStruct(myStruct,dirName,varargin)

p = inputParser;
addParameter(p,'useAlternative',false);
parse(p,varargin{:})
useAlternative = p.Results.useAlternative;

if nargin == 1 || isempty(dirName); dirName = 'myStructs'; end
if useAlternative; dirName = fullfile(dirName,'useAlternative');  end

subj = myStruct(1).subjID;
patientNum = subj(end-1:end);
lastwarn('')

for kk = 1:length(myStruct)
    if isfield(myStruct,'sourceLoc') && isfield(myStruct(kk).sourceLoc,'TimeSeries')
        myStruct(kk).sourceLoc = rmfield(myStruct(kk).sourceLoc,'TimeSeries');
        % warnOnce('Not removing time series.'); 
    end
end


save(sprintf('%s/myStruct%s',dirName,patientNum),'myStruct','-v7.3');

if ~isempty(lastwarn); warning('Apparent save failure for %s.',subj); failureFlag = true; 
else; fprintf('Save success for %s. \n',subj)
end

end