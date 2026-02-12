function traverseAllDirectories


directories = {'/Users/diamondjm/Documents/Josh_work/sourceLocSpikes/2022Jan26/myStructs';...
    '/Users/diamondjm/Documents/Josh_work/sourceLocSpikes/2022Jan26/myStructs/localizedStruct';...
    '/Users/diamondjm/Documents/Josh_work/sourceLocSpikes/2022Jan26/myStructs_depth';...
    '/Users/diamondjm/Documents/Josh_work/sourceLocSpikes/2022Jan26/myStructs_unified';...
    '/Users/diamondjm/Documents/Josh_work/sourceLocSpikes/2022Jan26/myStructs_unified/localizedStruct';...
    '/Users/diamondjm/Documents/Josh_work/sourceLocSeizure/2021Nov29/myStructs'};


for jj = 1:length(directories)
    patients = determinePatientsFromDirectory(directories{jj});

    for ii = 1:length(patients)
        
        [myStruct,failureFlag] = loadStruct(patients(ii),directories{jj});
        assert(~failureFlag);
        for kk = 1:length(myStruct)
            if isfield(myStruct(kk).sourceLoc,'braindata')
                myStruct(kk).sourceLoc = rmfield(myStruct(kk).sourceLoc,'braindata');
            end
        end
        saveStruct(myStruct,directories{jj})
    
    end
end

    