function checkTal(subj, rootEEGdir, isBatch, rerun)
% checkTal(subj, root, [ isBatch=0, rerun=0 ]) outputs metrics and summary figure files
%
% This should be run after localization is complete and outputs files that indicated that
% localization is complete to the master subject-completion checker

if nargin < 3 || isempty(isBatch), isBatch = 0; end
if nargin < 4 || isempty(rerun), rerun = 0; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script usage for Mike
if nargin < 1
    error('checkTal(subj,root)')
    subj = 'NIH055';
end
if nargin < 2 
    error('checkTal(subj,root)')
    dirs = evalin('base','dirs');
    rootEEGdir = dirs.zloc;
end

BATCH = 0; % save files without any prompt
RERUN = 0; % only used with batches, force re-run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Z

BATCH = isBatch || BATCH;
RERUN = rerun || RERUN;

if ~BATCH
    fprintf('----------------------------------------------------------------------------------------------------\n');
    fprintf('| Running localization metrics (check.m)\n');
    fprintf('----------------------------------------------------------------------------------------------------\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading + brain plot (f_metrics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bd              = braindata2();
bd.subj         = subj;
bd.rootEEGdir   = rootEEGdir;
bd.set_default_filepaths();

dmetrics = fullfile(bd.filepaths.zloc, 'metrics');
dtal = fileparts(bd.filepaths.zloc);
fname_summary = fullfile(dmetrics, 'metrics.txt');

if exist(fname_summary, 'file')
    
    if ~BATCH && ~inputYN(sprintf('%s found. Re-run check.m to regenerate metrics?', fname_summary))
        return
    end
    
    if BATCH && ~RERUN
        return
    end
end

if ~exist(dmetrics, 'dir')
    mkdir(dmetrics);
end

if ~exist(dtal, 'dir')
    mkdir(dtal);
end

bd.load_all_files();
bp = brainplotter();

fnames = fullfile(bd.filepaths.surfaces, {'std.141.lh.pial.gii','std.141.rh.pial.gii'});
bp.loadSurface(fnames, {'pial_lh', 'pial_rh'});

usedAlice = exist(fullfile(bd.filepaths.zloc, 'ALICE/3Dclustering'), 'dir') && ...
    numel(lsCell(fullfile(bd.filepaths.zloc, 'ALICE/3Dclustering'))) > 4;
leads = getLeads(subj, rootEEGdir);
proj = bd.tal.xyz;





f_metrics = figure('name','localization_summary', 'windowStyle','normal','units','normalized', 'outerposition',[0 0 1 1]);
subplot(2,2,1);
bp.plot({'pial_lh','pial_rh'});
bpax = gca;
bp.setOpacity(0.8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Anchor displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(bd.zloc, 'anchors')
    tanch = bd.zloc.anchors;
    d_anch = [];
    for i = 1:height(tanch)
        xyz0 = tanch{i,{'x','y','z'}};
        iproj = find(strcmpi(proj.chanName, tanch.chanName{i}), 1);
        xyz1 = proj{iproj,{'x','y','z'}};
        d = norm(xyz1 - xyz0);
        d_anch = [d_anch d];
        xyz = [xyz0; xyz1];
        line(xyz(:,1), xyz(:,2), xyz(:,3), 'lineWidth', 8, 'color','k');
        bp.plotPoint(xyz0, 'color', [1 1 1]);
    end
    isAnch = height(tanch) > 0;
else
    isAnch = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CT displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(bd.zloc, 'ct_xyz')
    ct = bd.zloc.ct_xyz;
    missing_ct = 0;
    utags = unique(util_split_stringnum(ct.chanName), 'stable');
else
    ct = table([],[],[],[],'VariableNames',{'chanName','x','y','z'});
    missing_ct = 1;
    tags = util_split_stringnum(proj.chanName);
    utags = unique(unique(tags));
end

ntags = numel(utags);
cols = hsv(ntags);
d_by_tag = cell(ntags,1);
d_all = [];
chans = [];
tags = [];
for i = 1:height(ct)
    itag = find(strcmpi(util_split_stringnum(ct.chanName{i}), utags), 1);
    c = mean([cols(itag,:); [1 1 1]]); % whitened color
    xyz0 = ct{i,{'x','y','z'}};
    bp.plotPoint(xyz0, 'color', c);
    iproj = find(strcmpi(proj.chanName, ct.chanName{i}), 1);
    xyz1 = proj{iproj,{'x','y','z'}};
    xyz = [xyz0; xyz1];
    line(xyz(:,1), xyz(:,2), xyz(:,3), 'lineWidth', 8, 'color','k');
    bp.plotPoint(xyz1, 'color', cols(itag,:), 'legend', utags{itag});

    chans = [chans ct.chanName(i)];
    tags = [tags utags(itag)];

    % distance
    d = norm(xyz1 - xyz0);
    d_by_tag{itag} = [d_by_tag{itag} d];
    d_all = [d_all d];
end

if missing_ct
    bp.plotPointGroup(proj, 'tagNames', tags);
end
suplabel(subj,'t');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boxplot of hardware CT displacement (with anchors) (f_subdurals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isAnch
    labels = [tags repmat({'Anchors'}, 1, height(tanch))];
    vals = [d_all d_anch];
else
    d_anch = 0;
    labels = tags;
    vals = d_all;
end
%f_subdurals = f(112);
subplot(2,2,2);
[a] = boxplot(vals, labels);
xlabel Hardware
set(gca,'XTickLabelRotation',45);
ylabel('Displacement (mm)');
title(subj);
figfmt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Histogram of CT displacement and dist2dura (f_multi_implant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dave fix 4.2019. Error here if patient only has depths. Adding if
% statement to avoid
if ~isempty(bd.zloc.dural_dist_1)

    d_dura = bd.zloc.dural_dist_1.dist2dural;
    %f_multi_implant = f(113); clf;
    subplot(2,2,3);
    histogram(d_all); hold on;
    histogram(d_dura);
    xlabel('Displacement (mm)');
    legend({'CT displacement' 'Dura distance'});
    figfmt
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Average spacing - CT and projected (f4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_nbr_ct    = [];
d_nbr_proj  = [];
bpchan = getLeads(subj,rootEEGdir, 'isBP',1);
chans = union(ct.chanName, proj.chanName);
for i = 1:length(bpchan)
    ichan_ct = find(strcmpi(ct.chanName, bpchan{i,1}), 1);
    jchan_ct = find(strcmpi(ct.chanName, bpchan{i,2}), 1);
    ichan_pr = find(strcmpi(proj.chanName, bpchan{i,1}), 1);
    jchan_pr = find(strcmpi(proj.chanName, bpchan{i,2}), 1);
    d_ct = norm(ct{ichan_ct,{'x','y','z'}} - ct{jchan_ct,{'x','y','z'}});
    d_proj = norm(proj{ichan_pr,{'x','y','z'}} - proj{jchan_pr,{'x','y','z'}});
    d_nbr_ct = [d_nbr_ct, d_ct];
    d_nbr_proj = [d_nbr_proj, d_proj];
end


if any(d_nbr_proj) >= 16.5
    fprintf('WARNING: inter-electrode spacing high!\n');
end

% neighbors by tag. Stored in matrix columns are each tag
% d_nbr_ct_hw = nan(length(bpchan), numel(utags));
% d_nbr_proj_hw = nan(length(bpchan), numel(utags));
hw = util_split_stringnum(bpchan);
% for i = 1:numel(utags)
%     temp = d_nbr_ct(strcmp(hw, utags{i}));
%     d_nbr_ct_hw(1:numel(temp), i) = temp;
%     
%     temp = d_nbr_proj(strcmp(hw, utags{i}));
%     d_nbr_proj_hw(1:numel(temp), i) = temp;
% end
subplot(2,2,4);
X = [d_nbr_ct, d_nbr_proj]';
G1 = [repmat({'CT'},size(d_nbr_ct)), repmat({'proj'},size(d_nbr_ct))]';
G2 = [hw; hw];
hierarchicalBoxplot(X, [categorical(G1) categorical(G2)]);
texts = findall(gca, 'type', 'text');
for i=1:numel(texts), texts(i).Parent = []; end
n = numel(unique(hw));
ax = gca;
tmp = ceil(ax.XLim(2) / n) * n;
ax.XTick = 1 : tmp/n : tmp;
ax.XTickLabel = unique(hw, 'stable');
ax.XTickLabelRotation = 45;

% remove boxplot's bottom label lines
i=1; 
while i <= numel(ax.Children)
    if strcmpi(ax.Children(i).Type, 'Line'), ax.Children(i).Parent = [];
    else, i=i+1;
    end
end
% reference line at typical spacing
hline(10);
hline(5);
%figfmt

%f4 = f(114); clf;

% histogram(d_nbr_ct); hold on;
% histogram(d_nbr_proj);
ylabel('Inter-electrode spacing (mm)');
%legend({'CT' 'Projected'});
figfmt

missingLeads = setdiff(leads, proj.chanName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multi-implant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(bd.zloc, 'ct_2_xyz')
    f_multi_implant = figure;
    bp.plot(fieldnames(bp.surfaces));
    bp.plotPoint(bd.zloc.proj_1_xyz, 'color',[0 1 0], 'legend','CT_1');
    bp.plotPoint(bd.zloc.proj_2_xyz, 'color',[0 0 1], 'legend','CT_2');
    bp.legend('orientation','horizontal');
    figfmt
    f_multi_implant.Name = 'mutli-implant';
    dual_ct = 1;
    title(subj);
else
    dual_ct = 0;
    f_multi_implant = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YN = {'NO', 'YES'};
yn = @(x) YN{x+1}; % 0/1 --> Y/N

fd = fopen(fname_summary, 'w+');
if fd > 0
    % settings
    fprintf(fd, '%-30s: %s\n', 'Used Anchors', yn(isAnch));
    fprintf(fd, '%-30s: %s\n', 'Used ALICE', yn(usedAlice));
    fprintf(fd, '%-30s: %s\n', 'Missing leads?', yn(numel(missingLeads)>0));
    for i = 1:numel(missingLeads)
        fprintf(fd, '%-23s\n', sprintf('\tMissing %s', missingLeads{i}));    
    end

    % spacing
    fprintf(fd, '%-30s: %.1f +/- %.1f (SEM)\n', 'Inter-electrode CT', mean(d_nbr_ct), std(d_nbr_ct)/sqrt(numel(d_nbr_ct)));
    fprintf(fd, '%-30s: %.1f +/- %.1f (SEM)\n', 'Inter-electrode Proj', mean(d_nbr_proj), std(d_nbr_proj)/sqrt(numel(d_nbr_proj)));

    % displacement
    fprintf(fd, '%-30s: %.1f +/- %.1f (SEM)\n', 'Anchor displacement', mean(d_anch), std(d_anch)/sqrt(numel(d_anch)));
    if missing_ct
        fprintf(fd, '%-30s: MISSING\n', 'CT');
    else
        
        fprintf(fd, '%-30s: %.1f +/- %.1f (SEM)\n', 'CT displacement', mean(d_all), std(d_all)/sqrt(numel(d_all)));
        for itag = 1:ntags
            d = d_by_tag{itag};
            fprintf(fd, '%-23s: %.1f +/- %.1f (SEM)\n', sprintf('\tCT diplacement %s',utags{itag}), mean(d), std(d)/sqrt(numel(d)));
        end
    end

    % dura
    % Dave fix 4.2019. Error here if patient only has depths. Adding if
    % statement to avoid
    if ~isempty(bd.zloc.dural_dist_1)
        
        fprintf(fd, '%-30s: %.1f +/- %.1f (SEM)\n', 'Distance from dura', mean(d_dura), std(d_dura)/sqrt(numel(d_dura)));
        fclose(fd);
        system(sprintf('cat %s', fname_summary));
        
    end
else
    warning('Could not write/open: %s', fname_summary);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Brain figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_subdurals = figure('Name', 'electrode_locations', 'windowStyle','normal', 'units','normalized','outerposition',[0 0 1 1]);
bp.viewAll(f_subdurals, bpax);
pause(5); % SJ
bp.legend('orientation', 'horizontal');
pause(5)
subl = suplabel(subj,'t');
pause(5)
set(subl,'FontSize',30);
pause(5); %SJ
figfmt(30,[]);
pause(5)


if BATCH || inputYN('Save figures?')
    figpdf(f_metrics,   [], dmetrics, 'png', 'overwrite',1 );
    figpdf(f_subdurals, [], dmetrics, 'png', 'overwrite',1);
    figpdf(f_subdurals, [], dtal, 'png', 'overwrite',1);
    if dual_ct
        figpdf(f_multi_implant, [], dmetrics, 'png', 'overwrite',1);
    end
    
    % Similar figure for depths
    depth_chans = getLeads(subj, rootEEGdir, 'hardwareType', {'depth'});
    f_depths = []; %#ok<NASGU>
    if ~isempty(depth_chans)
        t_depths = bd.tal.xyz(ismember(bd.tal.xyz.chanName, depth_chans), :);
        f_depths = figure('Name', 'electrode_locations_depth', 'windowStyle','normal', 'units','normalized','outerposition',[0 0 1 1]);
        axes(bpax); pause(7);
        bp.clearPoints();
        bp.setOpacity(0.25);
        bp.plotPointGroup(t_depths, 'element_info',bd.docs.element_info);
        bp.viewAll(f_depths, bpax); pause(5)
        bp.legend('orientation', 'horizontal'); pause(5)
        suplabel(subj,'t');
        pause(5)
        figfmt(30,[]);
        pause(5);
        figpdf(f_depths, [], dmetrics, 'png', 'overwrite',1 );
        figpdf(f_depths, [], dtal, 'png', 'overwrite',1);
    end
    
end

if BATCH
    close(f_depths, f_multi_implant, f_subdurals, f_metrics);
end

end % end checkTal

%% Begin private functions (from mike toolbox)

function figfmt(fontSize, lineWidth)
% figfmt(fontSize,lineWidth)
% Basic formating

    if nargin < 1 || isempty(fontSize), fontSize = 16; end
    if nargin < 2, lineWidth = 2; end
    f = gcf;
    f.Color = 'w';
    
    ax = findall(f,'type','axes');
    for i=1:length(ax)
        axes(ax(i));
        ax(i).FontSize = fontSize;
        %axis tight
    end
    
    lgd = findall(gcf,'type','legend');
    for i=1:length(lgd)
        lgd(i).FontSize = fontSize;
    end
    
    if ~isempty(lineWidth)
        figall('line', @(l) set(l,'LineWidth',lineWidth), gcf)
    end
    
    %box off
end

function figpdf(f, subd, d, ext, varargin)
% Saves figure f to ~/Desktop/figures/[date] or a given dir or subdir
% FIGPDF(f, subd, d, 'svg'/'pdf'/'fig') saves f to d/subd as an ext file
% FIGPDF(gcf,subd) will save a png to ~/figures/[date]/subd
%
%
% Inputs
%   f - figure (Default gcf)
%   subd - name of subdirectory (good if you're saving a group, Default none)
%   d - directory (Default ~/Desktop/figures)
%   ext - extension (Default png)
%
% Optional Input
%   saveAlpha - uses openGL; will NOT outpuut vector graphic. Workaround
%               for fact that eps doesn't support transparency. Consider
%               using svg
%   overwrite - force overwrite without input
%
% Details:
%   Name of file is either the figure name or (if no name) the currentAxis title
%   If you specify a subd, no more than MAX_SUBD_FIGS (=15) will be saved.
%   The oldest files will be moved to Trash!
%
% See also: print
    MAX_SUBD_FIGS = 50; % delete older figs
    
    % Parse input
    ip = inputParser;
    if nargin < 1 || isempty(f), f = gcf; end 
    if nargin < 2 || isempty(subd), subd = []; end
    if nargin < 3 || isempty(d), d = sprintf('~/figures/%s', datestr(datetime,'yy_mm_dd')); end
    if nargin < 4 || isempty(ext)
        ext='png'; fmt = '-dpng'; %'-dpdf';
        %figpdf(f, [], d, 'fig', varargin{:});
        %figpdf(f, [], d, 'png', varargin{:});
        %return;
    else
        fmt = ['-d' ext];
    end
    ip.addParameter('saveAlpha',0);
    ip.addParameter('overwrite',[]);
    ip.parse(varargin{:});
    opengl = ip.Results.saveAlpha;
    overwrite = ip.Results.overwrite;
    
    if numel(f) > 1
        for fi = 1:numel(f)
            ffi = f(fi);
            if isnumeric(ffi)
                ffi = figure(ffi); pause(0.5);
            end  
            figpdf(ffi, subd, d, ext, varargin{:}, 'overwrite',1);
        end
        return
    end
    
    is_subd = false;
    if ~isempty(subd) 
        d = fullfile(d,subd);
        is_subd = true;
        
    end
    
    % Makes sure the scale is good - this seems to produce a nice large 1:1 apect for giza
    f.PaperPositionMode = 'auto';
    f.Units = 'pixels';
    %f.Position = [1 214 1582 1131];
    
    
    % Create directory and check contents. Delete old figs if subd
    if ~exist(d,'dir'), mkdir(d); end
    figs = dir(d);
    figs = figs(cellfun(@(x) x(1)~='.', {figs.name}));
    [~,dateNdx] = sort([figs.datenum]);
    figs = figs(dateNdx);
    if is_subd
        fprintf('Warning: %d/%d figures (using a subd limits you to %d figures. Oldest figs will be deleted)\n',...
            numel(figs), MAX_SUBD_FIGS, MAX_SUBD_FIGS);
        if numel(figs) >= MAX_SUBD_FIGS
            delete(fullfile(d, figs(1).name));
        end
    end
    
    
    
    % Infer the filename from the figure name or title. If no name, use
    % an incremented number based on current figs in the directory
    name = f.Name;
    if isempty(name), try name = f.CurrentAxes.Title.String; catch, end, end
    if isempty(name)
    
        figs = lsCell(d);
        [names,nums] = util_split_stringnum(figs);
        nums = nums(strncmpi('fig',names,3));
        name = sprintf('fig%d',max([nums;0])+1);
    end
    filename = fullfile(d, name);
    
    % Check for existing file and overwrite if needed
    while exist(sprintf('%s.%s',filename,ext), 'file') && isempty(overwrite)
        if ~inputYN('file exists, overwrite?')
            filename = fullfile(d, input('Enter new name: ', 's'));
            overwrite = 0;
        else
            overwrite = 1;
        end
        if isempty(filename), return; end
    end
    
    
    % Call print or save
    if strcmp(ext, 'fig')
        savefig(f, filename);
    elseif opengl
        fprintf('Using openGL; output is *not* in vector format\n');
        print(f, filename ,fmt,'-r300', '-opengl');
    else
        try
            print(f, filename ,fmt,'-r0');%, '-bestfit');
        catch
            print(f, filename, fmt, '-r0');
        end
            
    end
    fprintf('Fig saved as: %s.%s\n', filename, ext);
end

function figall(type, fun, root, varargin)
% FIGALL(type, fun, groot, ...)
% iterates/selects each obj=TYPE and then calls the function fun(obj, varargin)
% example: figall('axes', @(obj) obj.CLim=[0,10])
    if nargin < 3, root = groot; end
    
    x = findall(root, 'type', type);
    for i = 1:length(x)
        
        if isprop(x(i), 'UserData') && isfield(x(i).UserData, 'skipFigall') && x(i).UserData.skipFigall
            continue;
        end
        
        if strcmpi(type,'axes')
            axes(x(i)); 
        end
        fun(x(i), varargin{:});
    end
end

function obj = hline(y, textLabel, varargin)
% plots a horizontal line at the given y coordinate. If textLabel=1, writes the value above line
% vline(y, textLabel, ...)
% plot([y y], [ax.XLim], varargin{:})
    hold on
    if nargin < 2, textLabel = 0; end
    ax = gca;
    obj = plot([ax.XLim], [y y], varargin{:});
    
    if textLabel
        text(ax.XLim(2), y, num2str(y));
    end
end