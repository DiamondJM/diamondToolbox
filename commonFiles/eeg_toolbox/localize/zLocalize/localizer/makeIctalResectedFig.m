function fig = makeIctalResectedFig(subj, rootEEGdir, type, rerun)
    % fig = makeIctalResectedFig(subj, rootEEGdir, [type]) makes a brain figure with ictal, interictal, and resected electrodes labeled
    %
    % INPUT:
    %   subj
    %   rootEEGdir
    %   type        - [OPTIONAL] One of ictal, interictal, or resected. If not passed, function is called once for each
    %   rerun       - [OPTIONAL] If 0 (or omitted) AND if figure already exists in expected file location, do nothing. 1 to rerun
    %
    %
    
    % OUTPUT:
    %   fig         - handle to figure of given type. If type is not given, no handles are output.
    %
    % FILE OUTPUT:
    %   figs are saved to subject docs directory as png's
    %
    % REVISION HISTORY:
    %   04/18 MST   - Created
    
    fig = [];
    
    if ~exist(fullfile(rootEEGdir,subj,'tal/leads.csv'),'file'),
        fprintf('\n HEADS UP: tal/leads.csv folder not created, so cant make ictal/interictal/resected figures yet');
        return;
    end
    
    if nargin < 4 || isempty(rerun)
        rerun = 0;
    end

    
    if nargin < 3 || isempty(type)
        makeIctalResectedFig(subj, rootEEGdir, 'ictal', rerun);
        makeIctalResectedFig(subj, rootEEGdir, 'interictal', rerun);
        makeIctalResectedFig(subj, rootEEGdir, 'resected', rerun);
        return
    else
        fig = handle([]);
    end
    
    
    SETTINGS.BRAIN_OPACITY = 0.6;
    SETTINGS.BRAIN_OPACITY_WITH_DEPTHS = 0.35;
    SETTINGS.MARKED_OPACITY = 0.85;
    SETTINGS.NORMAL_OPACITY = 0.85;
    SETTINGS.OPACITY_DEPTH  = 1.0;
    COLORS.SURFACE  = [1 1 1] * 1.0;
    COLORS.SURFACE_MIX = 0.5;       % 1.0 gives full colors.surface value; 0 gives default gray brain
    COLORS.RESECTED = [1 0 0] * 1.0; 
    COLORS.ICTAL    = [1 0 0] * 1.0;
    COLORS.INTERICTAL = [1 0 0] * 1.0;
    COLORS.NORMAL   = [1 1 1] * 1.0;
    COLORS.TEXT     = [0 0 0];
    RADIUS          = 2.5;
    RADIUS_DEPTH    = 1.5;
    FILES.DIR       = fullfile(rootEEGdir, subj, 'docs');
    FILES.RESECTED  = fullfile(FILES.DIR, 'electrodes_resected.png');
    FILES.ICTAL     = fullfile(FILES.DIR, 'electrodes_ictal.png');
    FILES.INTERICTAL = fullfile(FILES.DIR, 'electrodes_interictal.png');
    
    % parse type
    % e.g. isResected --> remove the "is" 
    type = strrep(lower(type), 'is', '');
    assert(ismember(type, {'ictal', 'interictal', 'resected'}), 'invalid type %s', type);
    
    
    % check file
    fpath = FILES.(upper(type));
    [~,fname] = fileparts(fpath);
    if exist(fpath, 'file') && ~rerun
        fprintf('File already exists: %s\n', fpath);
        return
    end
    
    
    % parse the patient info file
    %   this will provide a written version of the resected tissue and indicates the rare cases when no ictal or interictal electrodes were recorded
    try
        patient_info = parsePatientInfo(subj, rootEEGdir, 0);
    catch e
        warning(e.identifier, '%s', e.message);
        fprintf('\n if this error still happens, then need to fix patient info, or come up with defaults for patient_info parsing');
        patient_info.resection = sprintf('<<Unable to parse patient info>>');
        keyboard;
    end

    noneStr = '';
    chanNames = getLeads(subj, rootEEGdir, 'markedAs', type);
    if isempty(chanNames)
        expectNoChan = 0;
        
        if strcmp(type,'resected'),
            if contains(lower(patient_info.resection),'none'),
               expectNoChan = 1;
            end
        elseif strcmp(type,'ictal') & patient_info.Ictal_isNone==1,
            expectNoChan = 1;
            noneStr = 'None Found During Monitoring';
        elseif strcmp(type,'interictal') & patient_info.Interictal_isNone==1,
            expectNoChan = 1;
            noneStr = 'None Found During Monitoring';
        end
        
        %- channels not identified in element_info, and patient info does not note as "none", so there is a discrepancy
        if expectNoChan==0,
            fprintf('\n\nWARNING:  no channels marked %s in element info, but patient info doesnt indicate "%s:\r none".  \n Possibly kareem just hasnt updated these files yet.  \n Fix one or the other (or email kareem) and rerun makeIctalResectedFig',type,type);
            %- if a "none" is missing from patient_info, then add the break, "Ictal:" and "None" "------------------\nIctal:\nNone"
            keyboard;
            return
        end
    end
    
    
    % Plot surface
    bd = braindata2(subj, rootEEGdir);
    bp = bd.ezplot([],[],1); %- 12/2018, JW added forth (optional) param makes explot pick the first of multiple options if ambiguous
    bp.setOpacity(SETTINGS.BRAIN_OPACITY);
    bp.clearPoints();
    pause(2);
    fig_base = gcf();
    
    surfNames = fieldnames(bp.surfaces);
    for i=1:numel(surfNames)
        pause(0.5) %- occasional crash on the plotRegion step... see if this helps
        surfName = surfNames{i};
        bp.plotRegion(1:bp.stdNumVert, 'color', COLORS.SURFACE, 'surf',surfName, 'blend',COLORS.SURFACE_MIX);
    end
    
    
    col = COLORS.(upper(type));
    
    
    jack = bd.docs.jacksheet;
    subdurals = jack.chanName(strcmpi(jack.hardwareType, 'SUBDURAL'));
    assert(isprop(bd, 'tal') && isfield(bd.tal, 'xyz') && ~isempty(bd.tal.xyz), 'No bd.tal.xyz data. Possible cause: localization not run yet');
    xyz = bd.tal.xyz(ismember(bd.tal.xyz.chanName, intersect(subdurals, chanNames)),:);
    [~,nums] = util_split_stringnum(xyz.chanName);
    label = arrayfun(@num2str, nums,'uniformOutput',0);
    bp.plotPoint(xyz, 'color', col, 'alpha',SETTINGS.NORMAL_OPACITY, 'label',label, 'labelColor', COLORS.TEXT, 'radius',RADIUS);
    
    xyz_other =  bd.tal.xyz(ismember(bd.tal.xyz.chanName, setdiff(subdurals, chanNames)),:);
    bp.plotPoint(xyz_other, 'color', COLORS.NORMAL, 'alpha',SETTINGS.NORMAL_OPACITY, 'radius',RADIUS);
    
    %- add depth plotting
    depths    = jack.chanName(strcmpi(jack.hardwareType, 'DEPTH'));
    xyz_depth = bd.tal.xyz(ismember(bd.tal.xyz.chanName, intersect(depths, chanNames)),:);
    bp.plotPoint(xyz_depth, 'color', col, 'alpha',SETTINGS.OPACITY_DEPTH, 'radius',RADIUS_DEPTH);
    
    xyz_depth_other = bd.tal.xyz(ismember(bd.tal.xyz.chanName, setdiff(depths, chanNames)),:);
    bp.plotPoint(xyz_depth_other, 'color', COLORS.NORMAL, 'alpha',SETTINGS.OPACITY_DEPTH, 'radius',RADIUS_DEPTH);
    
    if length(depths)>0 && length(intersect(depths, chanNames)),
        bp.setOpacity(SETTINGS.BRAIN_OPACITY_WITH_DEPTHS);
    end
    
    pause(5);
    
    
    
    
    % Figure is how we want it. Copy it to a new figure in multiple views.
    bpax = gca();
    fig = figure('name', fname, 'windowStyle','normal','units','normalized', 'outerposition',[0 0 1 1]); 
    pause(1);
    axes(fig.CurrentAxes); 
    pause(1);
    bp.viewAll(fig, bpax);
    pause(5);
    
    
    
    % Add a subtitle for resected figure
    numChans = numel(getLeads(subj, rootEEGdir));
    if strcmpi(type, 'resected')
        titleStrExtra = sprintf('%d of %d electrodes resected\n%s', numel(chanNames), numChans, patient_info.resection);
        
        if numel(chanNames)==0 & ~contains(lower(patient_info.resection),'none'),
            fprintf('\n ERROR:  resection specified in patient info but no electrodes included in element_info. Must update element_info!!!!\NOT SAVING FIG!\n');
            %if ~BATCH, keyboard; end
            close(fig);
            close(fig_base);
            delete(fpath);  %- remove the fig that is already there if previously saved
            return;
        end
        
    else
        titleStrExtra = sprintf('%d of %d electrodes identified\n%s', numel(chanNames), numChans, noneStr);
    end
    
    if isempty(titleStrExtra)
        title_str = sprintf('%s %s', subj, type);
    else
        title_str = sprintf('%s %s\n%s', subj, type, titleStrExtra);
    end
    
    %- fig title and font/etc settings
    [tax,txt] = suplabel(title_str, 't');
    figfmt(25,[]);
    
    
    % Save to file
    figpdf(fig, [], fileparts(fpath), 'png', 'overwrite',1);
    pause(2)
    
    close(fig);
    close(fig_base);
end



%% Begin private functions (from mike toolbox)
function figfmt(fontSize, lineWidth)
% figfmt(fontSize,lineWidth)
% Basic formating

    if nargin < 1 || isempty(fontSize), fontSize = 15; end
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
        %lgd(i).FontSize = fontSize;
        lgd(i).Visible = 'off';  %- for resected/interictal figs just turn this off
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