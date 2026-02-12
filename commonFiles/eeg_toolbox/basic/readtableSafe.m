function table = readtableSafe(filename, varargin)
% READTABLESAFE calls readtable with ','-delimiter and catches common errors
%
% Error-handling: 
%   UnequalVarLengths - adds a newline character to the end of file to fix
%
% Argument changes:
%   MATLAB2016a has a problem detecting ','-delimiter, so use this is delim
%   param not passed
%
% 6/10/2020 SJ: MATLAB 2020a does not seem to read some tables in correctly (like element_info_key), so
%               just run with 'Format', 'auto', which runs it with the 2019b version
% 7/7/2020 SJ: If you are trying to read in a spreadsheet, just run normally without adding delimiter or format option
% 

    newArgs = {'delimiter', ','};
    

    warning('off','MATLAB:textscan:UnableToGuessFormat');
    warning('off','MATLAB:textscan:AllNatSuggestFormat');
    mat_ver = version('-release');
    
    if contains(class(detectImportOptions(filename)),'Spreadsheet') %XLS or XLSX files
        % This is a spreadsheet that does not accept the 'format' or 'delimiter' parameters, do not want to add newArgs
        table = readtable(filename,varargin{:});
    else %CSV or TXT files
        varargin = [newArgs, varargin]; % note parser takes last arg, so you can override this if you want
        try
            if str2double(mat_ver(1:4))>2019 %SJ: Fix inconsistencies with 2020 version of readtable
                table = readtable(filename, 'Format','auto', varargin{:});
            else
                table = readtable(filename, varargin{:});
            end
        catch ME
            if strcmp(ME.identifier, 'MATLAB:readtable:UnequalVarLengthsFromFileNoFormat') || ...
                    strcmp(ME.identifier, 'MATLAB:readtable:UnequalVarLengthsFromFileWithFormat')
                % for some reason, adding a newline often fixes a table read error
                system(sprintf('echo "" >> %s', filename));
                table = readtable(filename, varargin{:});

            elseif strcmp(ME.identifier, 'MATLAB:FileIO:InvalidFid')
                % Someone is forgetting to fclose somewhere
                fclose all;
                table = readtable(filename, varargin{:});
            else
                rethrow(ME);
            end
        end
    end
    warning('on','MATLAB:textscan:UnableToGuessFormat');
    warning('on','MATLAB:textscan:AllNatSuggestFormat');
end