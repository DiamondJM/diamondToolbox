function sdAnnotations = readSdAnnotations(subj, rootFolder, varargin)
% readSdAnnotations  Parse SD annotation CSV for a subject.
%
%   sdAnnotations = readSdAnnotations(subj, rootFolder)
%   sdAnnotations = readSdAnnotations(subj, rootFolder, 'chanNames', cNames)
%
%   Returns a struct array. Each element is one SD group (one 'SD' row):
%     .groupTime   datetime     onset of the SD group
%     .electrodes  struct array with fields:
%                    .name   char      electrode label (from chanNames if given,
%                                      otherwise the raw CSV label e.g. 'SD6')
%                    .time   datetime  annotated onset for this electrode

ip = inputParser;
ip.addRequired('subj',       @ischar);
ip.addRequired('rootFolder', @ischar);
ip.addParameter('chanNames', {}, @(x) iscell(x) || isstring(x));
ip.parse(subj, rootFolder, varargin{:});
chanNames = cellstr(ip.Results.chanNames);
if isequal(chanNames, {''})
    chanNames = {};
end

csvPath = fullfile(rootFolder, subj, 'ts', 'Annotations_SD.csv');
assert(isfile(csvPath), '[readSdAnnotations] File not found: %s', csvPath);

fid = fopen(csvPath, 'r');
C   = textscan(fid, '%q%q%*[^\n]', 'Delimiter', ',', 'HeaderLines', 10);
fclose(fid);
dtStrs = C{1};
names  = C{2};

dtFmt = 'M/d/yyyy H:mm';
dts   = NaT(numel(dtStrs), 1);
for ii = 1:numel(dtStrs)
    try
        dts(ii) = datetime(dtStrs{ii}, 'InputFormat', dtFmt);
    catch
    end
end

sdAnnotations = struct('groupTime', cell(0,1), 'electrodes', cell(0,1));
currentGroup  = [];

for ii = 1:numel(names)
    n = strtrim(names{ii});
    if strcmp(n, 'SD')
        if ~isempty(currentGroup)
            sdAnnotations(end+1, 1) = currentGroup; %#ok<AGROW>
        end
        currentGroup = struct('groupTime', dts(ii), ...
            'electrodes', struct('name', {}, 'time', {}));
    elseif ~isempty(currentGroup) && ~isnat(dts(ii))
        tok = regexp(n, '^SD(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            elecNum  = str2double(tok{1});
            elecName = resolveElecName(n, elecNum, chanNames);
            currentGroup.electrodes(end+1) = struct('name', elecName, 'time', dts(ii));
        end
    end
end
if ~isempty(currentGroup)
    sdAnnotations(end+1, 1) = currentGroup;
end

end

% -------------------------------------------------------------------------
function name = resolveElecName(rawLabel, elecNum, chanNames)
if isempty(chanNames)
    name = rawLabel;
    return
end
pat = sprintf('(?<!\\d)%d(?!\\d)', elecNum);
for ii = 1:numel(chanNames)
    if ~isempty(regexp(chanNames{ii}, pat, 'once'))
        name = chanNames{ii};
        return
    end
end
name = rawLabel;
end
