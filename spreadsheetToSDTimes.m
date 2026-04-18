function [sdAnnotations, clipStart] = spreadsheetToSDTimes(rootFolder, subj, varargin)
% spreadsheetToSDTimes  Parse SD annotation CSV for a subject.
%
%   sdAnnotations = spreadsheetToSDTimes(rootFolder, subj)
%   sdAnnotations = spreadsheetToSDTimes(rootFolder, subj, 'chanNames', cNames)
%
%   Returns a struct array. Each element is one SD group (one 'SD' row):
%     .groupTime     datetime   absolute clock time of the SD group marker
%     .groupClipTime duration   time from clip start to group marker
%     .durationOverall    double     SD group duration in seconds (column C)
%     .electrodes    struct array with fields:
%                      .name      char      electrode label
%                      .time      datetime  absolute clock time
%                      .clipTime  duration  time from clip start to detection

ip = inputParser;
ip.addRequired('rootFolder', @ischar);
ip.addRequired('subj',       @ischar);
ip.addParameter('chanNames', {}, @(x) iscell(x) || isstring(x));
ip.parse(rootFolder, subj, varargin{:});
chanNames = cellstr(ip.Results.chanNames);
if isequal(chanNames, {''})
    chanNames = {};
end

tsFolder = fullfile(rootFolder, subj, 'ts');
hits = dir(fullfile(tsFolder, '*.csv'));
hits = hits(~cellfun(@isempty, regexpi({hits.name}, 'annotations', 'once')));
assert(~isempty(hits), '[spreadsheetToSDTimes] No annotations CSV found in %s', tsFolder);
csvPath = fullfile(tsFolder, hits(1).name);

% Dynamically find where data rows start (first line beginning with a date M/D/...)
fid = fopen(csvPath, 'r');
headerLines = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && ~isempty(regexp(line, '^\d{1,2}/\d{1,2}/', 'once'))
        break
    end
    headerLines = headerLines + 1;
end
fclose(fid);

fid = fopen(csvPath, 'r');
C   = textscan(fid, '%q%q%q%*[^\n]', 'Delimiter', ',', 'HeaderLines', headerLines);
fclose(fid);
dtStrs    = C{1};
names     = C{2};
durations = str2double(C{3});

dtFmt = 'M/d/yyyy H:mm';
dts   = NaT(numel(dtStrs), 1);
for ii = 1:numel(dtStrs)
    try
        dts(ii) = datetime(dtStrs{ii}, 'InputFormat', dtFmt);
    catch
    end
end

% Find clip start from first 'Start Recording' annotation
clipStart = NaT;
for ii = 1:numel(names)
    if strcmp(strtrim(names{ii}), 'Start Recording')
        clipStart = dts(ii);
        break
    end
end
assert(~isnat(clipStart), '[spreadsheetToSDTimes] No "Start Recording" row found in %s', csvPath);

sdAnnotations = struct('groupTime',{}, 'groupClipTime',{}, 'durationOverall',{}, 'electrodes',{});
currentGroup  = [];

for ii = 1:numel(names)
    n = strtrim(names{ii});
    if strcmp(n, 'SD')
        if ~isempty(currentGroup)
            sdAnnotations(end+1, 1) = currentGroup; %#ok<AGROW>
        end
        currentGroup = struct( ...
            'groupTime',     dts(ii), ...
            'groupClipTime', dts(ii) - clipStart, ...
            'durationOverall',    seconds(durations(ii)), ...
            'electrodes',    struct('name',{},'time',{},'clipTime',{}));
    elseif ~isempty(currentGroup) && ~isnat(dts(ii))
        tok = regexp(n, '^SD(\d+)$', 'tokens', 'once');
        if ~isempty(tok)
            elecNum  = str2double(tok{1});
            elecName = resolveElecName(n, elecNum, chanNames);
            currentGroup.electrodes(end+1) = struct( ...
                'name',     elecName, ...
                'time',     dts(ii), ...
                'clipTime', dts(ii) - clipStart);
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
