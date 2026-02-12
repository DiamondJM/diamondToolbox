function fig2pngSimple(figHandle,figFile)
%
% Simple export fig to png
%  figHandle  - pass the value returned by gcf (get current figure)
%  figFile    - destination path for output png
%
%  8/2017 JW pulled this from "export_fig"... which is complicated and has lots of options and files
%       guts below come from exoprt_fig.m and print2array.m
%
 

%%%% BUNCH OF STUFF TO GET IMAGE READY

% Make sure the figure is rendered correctly _now_ so that properties like
% axes limits are up-to-date.
drawnow;
pause(0.05);  % this solves timing issues with Java Swing's EDT (http://undocumentedmatlab.com/blog/solving-a-matlab-hang-problem)

% Parse the input arguments
if isgraphics(figHandle),
    fig = figHandle;
else
    fig = get(0, 'CurrentFigure');
end
%[fig, options] = parse_args(nargout, fig, varargin{:});
options.magnify = 1;  %options.aa_factor=1; %- JW

% Ensure that we have a figure handle
if isequal(fig,-1)
    return;  % silent bail-out
elseif isempty(fig)
    error('No figure found');
end

% Isolate the subplot, if it is one
cls = all(ismember(get(fig, 'Type'), {'axes', 'uipanel'}));
if cls
    % Given handles of one or more axes, so isolate them from the rest
    fig = isolate_axes(fig);
else
    % Check we have a figure
    if ~isequal(get(fig, 'Type'), 'figure');
        error('Handle must be that of a figure, axes or uipanel');
    end
    % Get the old InvertHardcopy mode
    old_invhardcopy = get(fig, 'InvertHardcopy');
end

% % Hack the font units where necessary (due to a font rendering bug in print?).
% % This may not work perfectly in all cases.
% % Also it can change the figure layout if reverted, so use a copy.
% magnify = options.magnify * options.aa_factor;
% if isbitmap(options) && magnify ~= 1
%     fontu = findall(fig, 'FontUnits', 'normalized');
%     if ~isempty(fontu)
%         % Some normalized font units found
%         if ~cls
%             fig = copyfig(fig);
%             set(fig, 'Visible', 'off');
%             fontu = findall(fig, 'FontUnits', 'normalized');
%             cls = true;
%         end
%         set(fontu, 'FontUnits', 'points');
%     end
% end

try
    % MATLAB "feature": axes limits and tick marks can change when printing
    Hlims = findall(fig, 'Type', 'axes');
    if ~cls
        % Record the old axes limit and tick modes
        Xlims = make_cell(get(Hlims, 'XLimMode'));
        Ylims = make_cell(get(Hlims, 'YLimMode'));
        Zlims = make_cell(get(Hlims, 'ZLimMode'));
        Xtick = make_cell(get(Hlims, 'XTickMode'));
        Ytick = make_cell(get(Hlims, 'YTickMode'));
        Ztick = make_cell(get(Hlims, 'ZTickMode'));
    end
    
    % Set all axes limit and tick modes to manual, so the limits and ticks can't change
    % Fix Matlab R2014b bug (issue #34): plot markers are not displayed when ZLimMode='manual'
    set(Hlims, 'XLimMode', 'manual', 'YLimMode', 'manual');
    set_tick_mode(Hlims, 'X');
    set_tick_mode(Hlims, 'Y');
    if ~using_hg2(fig)
        set(Hlims,'ZLimMode', 'manual');
        set_tick_mode(Hlims, 'Z');
    end
catch
    % ignore - fix issue #4 (using HG2 on R2014a and earlier)
end

% Fix issue #21 (bold TeX axes labels/titles in R2014b)
try
    if using_hg2(fig)
        % Set the FontWeight of axes labels/titles to 'normal'
        texLabels = findall(fig, 'type','text', 'FontWeight','bold');
        set(texLabels, 'FontWeight','normal');
    end
catch
    % ignore
end

% Fix issue #42: non-normalized annotations on HG1 (internal Matlab bug)
annotationHandles = [];
try
    if ~using_hg2(fig)
        annotationHandles = findall(fig,'Type','hggroup','-and','-property','Units','-and','-not','Units','norm');
        originalUnits = get(annotationHandles,'Units');
        set(annotationHandles,'Units','norm');
    end
catch
    % should never happen, but ignore in any case - issue #50
end

% Fix issue #46: Ghostscript crash if figure units <> pixels
oldFigUnits = get(fig,'Units');
set(fig,'Units','pixels');

% Set to print exactly what is there
set(fig, 'InvertHardcopy', 'off');



%%%%%%%%

% now " [A, tcol] = print2array(fig, magnify, renderer);"   PEEL IT APART HERE:
res = 1;


% Warn if output is large
old_mode = get(fig, 'Units');
set(fig, 'Units', 'pixels');
px = get(fig, 'Position');
set(fig, 'Units', old_mode);
npx = prod(px(3:4)*res)/1e6;
if npx > 30
    % 30M pixels or larger!
    warning('MATLAB:LargeImage', 'print2array generating a %.1fM pixel image. This could be slow and might also cause memory problems.', npx);
end
% Retrieve the background colour
bcol = get(fig, 'Color');
% Set the resolution parameter
res_str = ['-r' num2str(ceil(get(0, 'ScreenPixelsPerInch')*res))];
% Generate temporary file name
tmp_nam = [tempname '.tif'];
try
    % Ensure that the temp dir is writable (Javier Paredes 26/2/15)
    fid = fopen(tmp_nam,'w');
    fwrite(fid,1);
    fclose(fid);
    delete(tmp_nam);  % cleanup
    isTempDirOk = true;
catch
    % Temp dir is not writable, so use the current folder
    [dummy,fname,fext] = fileparts(tmp_nam); %#ok<ASGLU>
    fpath = pwd;
    tmp_nam = fullfile(fpath,[fname fext]);
    isTempDirOk = false;
end


if nargin < 3
    renderer = '-opengl';
end
err = false;
% Set paper size
old_pos_mode = get(fig, 'PaperPositionMode');
old_orientation = get(fig, 'PaperOrientation');
set(fig, 'PaperPositionMode', 'auto', 'PaperOrientation', 'portrait');
try
    % Print to tiff file
    print(fig, renderer, res_str, '-dtiff', tmp_nam);
    % Read in the printed file
    A = imread(tmp_nam);
    % Delete the temporary file
    delete(tmp_nam);
catch ex
    err = true;
end

% Reset paper size
set(fig, 'PaperPositionMode', old_pos_mode, 'PaperOrientation', old_orientation);
% Throw any error that occurred
if err
    % Display suggested workarounds to internal print() error (issue #16)
    fprintf(2, 'An error occured with Matlab''s builtin print function.\nTry setting the figure Renderer to ''painters'' or use opengl(''software'').\n\n');
    rethrow(ex);
end
% Set the background color
if isequal(bcol, 'none')
    bcol = [];
else
    bcol = bcol * 255;
    if isequal(bcol, round(bcol))
        bcol = uint8(bcol);
    else
        bcol = squeeze(A(1,1,:));
    end
end

% Check the output size is correct
if isequal(res, round(res))
    px = round([px([4 3])*res 3]);  % round() to avoid an indexing warning below
    if ~isequal(size(A), px)
        % Correct the output size
        A = A(1:min(end,px(1)),1:min(end,px(2)),:);
    end
end




%%%%%% BACK TO MAIN EXPORT FIG

res = options.magnify * get(0, 'ScreenPixelsPerInch') / 25.4e-3;

%imwrite(A, [options.name '.png'], 'ResolutionUnit', 'meter', 'XResolution', res, 'YResolution', res);
imwrite(A, figFile, 'ResolutionUnit', 'meter', 'XResolution', res, 'YResolution', res);


%%%%%
%%%%% And reset back to original figure settings
% Reset the hardcopy mode
set(fig, 'InvertHardcopy', old_invhardcopy);
% Reset the axes limit and tick modes
for a = 1:numel(Hlims)
    try
        set(Hlims(a), 'XLimMode', Xlims{a}, 'YLimMode', Ylims{a}, 'ZLimMode', Zlims{a}, 'XTickMode', Xtick{a}, 'YTickMode', Ytick{a}, 'ZTickMode', Ztick{a});
    catch
        % ignore - fix issue #4 (using HG2 on R2014a and earlier)
    end
end
% Revert the tex-labels font weights
try set(texLabels, 'FontWeight','bold'); catch, end
% Revert annotation units
for handleIdx = 1 : numel(annotationHandles)
    try
        oldUnits = originalUnits{handleIdx};
    catch
        oldUnits = originalUnits;
    end
    try set(annotationHandles(handleIdx),'Units',oldUnits); catch, end
end
% Revert figure units
set(fig,'Units',oldFigUnits);



