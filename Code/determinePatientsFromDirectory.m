function allPatients = determinePatientsFromDirectory(chosenDirectory,useStructs)

if nargin == 1; useStructs = true; end

if useStructs
    if contains(chosenDirectory,'myStructs'); fn = chosenDirectory; 
    else; fn = fullfile(chosenDirectory,'myStructs');
    end
else
    fn = chosenDirectory;
end
fromDirectory = lsCell(fn);

for ii = 1:length(fromDirectory)
    fromDirectory{ii} = regexp(fromDirectory{ii},'\d*','Match');
end

fromDirectory(cellfun('isempty',fromDirectory)) = []; 
allPatients = zeros(size(fromDirectory));

for ii = 1:length(fromDirectory)
    allPatients(ii) = str2double(fromDirectory{ii});
end

allPatients = unique(allPatients);

end