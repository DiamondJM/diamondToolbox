function startTime = readBr(subj,sess)

% Intended for BR

subjInfo = getMicroSubjInfo_v11(subj);
rawPath = fullfile(subjInfo.dataPath56,subj,'data_raw',sess,'/LogFile_ieeg.txt'); 
if ~exist(rawPath,'file')
    rawPaths = lsCell(fullfile(subjInfo.dataPath56,subj,'data_raw',sess)); 
    rawPathsLog = contains(rawPaths,'LogFile'); 
    if sum(rawPathsLog)==0
        % Finally, just read in filename? 
        fprintf('Couldn''t find any log files, attempting to read in from filename. \n'); 
        
        [~,fn] = fileparts(rawPaths); 
        fn = fn(contains(fn,'INST'));
        if isempty(fn)
            [~,fn] = fileparts(rawPaths);
            fn = fn(contains(fn,{'INST','ieeg','utah'}));
            assert(isequal(sess(1:11),fn{1}(1:11)))
            fprintf('Couldn''t read seconds from filename, will simply pass back session. \n');
            fprintf('Should double check that seizures start at the right time. \n'); 
            startTime = datetime(sess,'InputFormat','yyMMdd_HHmm'); 
            return
            
        else
            
            fn = cellfun(@(x) strsplit(x,'-'),fn,'UniformOutput',0);
            fnDate = cellfun(@(x) x{2},fn,'UniformOutput',0);
            fnTime = cellfun(@(x) x{3},fn,'UniformOutput',0);
            % fnTime = cellfun(@(x) str2double(x), fnTime);
            
            assert(length(unique(fnDate)) == 1);
            assert(isequal(sess(1:6),fnDate{1}(3:end)));
            % Date is good, let's check time
            
            sessDate = convertFoldersToDates(sess);
            fnTimeDt = NaT(size(fnTime));
            for ii = 1:length(fnTime)
                fnTimeDt(ii) = datetime([sess(1:6) fnTime{ii}],'InputFormat','yyMMddHHmmss');
            end
            assert(all(fnTimeDt-sessDate < minutes(1)));
            if length(unique(fnTime)) > 1
                assert(range(fnTimeDt) < seconds(2));
                warning('Multiple second start times found in filenames. Choosing latest.');
            end
            fnTimeDt = max(fnTimeDt);
            
            startTime = fnTimeDt;
            return
        end
        
    else
        if sum(rawPathsLog)==2
            rawPathsLog = find(rawPathsLog);
            rawPathsLog = rawPathsLog(1);
            warning('Multiple log files found. Choosing %s arbitrarily.',rawPaths{rawPathsLog});
        end
        rawPath = fullfile(subjInfo.dataPath56,subj,'data_raw',sess,rawPaths{rawPathsLog});
    end
end

opts = detectImportOptions(rawPath,'NumHeaderLines',0);
opts.VariableNames = {'Date','Time','Notes'};
opts = setvaropts(opts,'Date','InputFormat','MM/dd/uuuu');
opts.Delimiter = {'\t'};

t = readtable(rawPath,opts);
assert(isequal(t.Notes{1},'File: Patient Data Save started Programmatically!'));

sessDate = datetime(sess,'InputFormat','yyMMdd_HHmm'); 

[y,m,d] = ymd(t.Date(1)); [yy,mm,dd] = ymd(sessDate);
assert(isequal([y m d],[yy mm dd]));
[h,m,s] = hms(t.Time(1)); [hh,mm] = hms(sessDate);
% assert(isequal([h m],[hh mm]));
assert(hours(h) + minutes(m) + seconds(s) - hours(hh) - minutes(mm) < minutes(5))

% startTime = sessDate + seconds(s);
startTime = t.Time(1) + t.Date(1);
startTime.Format = 'dd-MMM-uuuu HH:mm:ss';


end