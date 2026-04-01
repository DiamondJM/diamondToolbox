function ctCoregNii = ctMrRegistrationGUI(mrNii, ctNii, outputDir)
% CTMRREGISTRATIONGUI  Interactive rigid-body CT-MR registration GUI.
%
% Usage:
%   ctCoregNii = ctMrRegistrationGUI(mrNii, ctNii, outputDir)
%
% Inputs:
%   mrNii      - path to MRI NIfTI (.nii)
%   ctNii      - path to CT NIfTI (.nii)
%   outputDir  - directory to write ct_manual_coreg.nii
%
% Output:
%   ctCoregNii - path to saved file, or '' if cancelled
%
% Interaction:
%   Mode buttons: Translate | Rotate | Scale
%   Click+drag in any of the 3 views to apply that mode's transform.
%     Translate  – axial drag moves x/y; coronal x/z; sagittal y/z
%     Rotate     – left/right drag rotates around view normal (z/y/x)
%     Scale      – left/right drag anywhere scales isotropically
%   Click without drag: move crosshair
%   Scroll: advance through slices in the view under cursor
%   Shift: fine-adjust mode (1/10 sensitivity)
%   Save: resamples CT onto MRI grid, writes ct_manual_coreg.nii
%   Cancel / close: returns ''

%% ---- Load volumes ----
mrInfo = niftiinfo(mrNii);
mrVol  = double(niftiread(mrInfo));
ctInfo = niftiinfo(ctNii);
ctVol  = double(niftiread(ctInfo));

[nxMr, nyMr, nzMr] = size(mrVol);
pdMr = mrInfo.PixelDimensions(1:3);   % [dx dy dz] mm/vox

% NIfTI affines, row-vector convention: [vox_1based, 1] * T = [world, 1]
Tmr = mrInfo.Transform.T;
Tct = ctInfo.Transform.T;

%% ---- State ----
tx = 0; ty = 0; tz = 0;    % translation mm
rx = 0; ry = 0; rz = 0;    % rotation deg
sc = 1.0;                   % isotropic scale
ctAlpha = 0.5;              % CT overlay opacity [0-1]
curVox = [round(nxMr/2), round(nyMr/2), round(nzMr/2)];
mode_ = 'translate';        % 'translate' | 'rotate' | 'scale'

mrWW = 1000; mrWL = 500;    % MRI window width / level (brain)
ctWW = 3000; ctWL = 700;    % CT  window width / level (bone)

% CT centre in world (rotation/scale pivot)
ctCtrVox = [(size(ctVol,1)+1)/2, (size(ctVol,2)+1)/2, (size(ctVol,3)+1)/2];
ctCtrWorld = [ctCtrVox, 1] * Tct;   % 1×4

% Drag state
isDragging   = false;
dragFigStart = [];
dragState0   = [];          % [tx ty tz rx ry rz sc] at drag-down
dragAxSel    = [];          % which axes was clicked

% Sensitivity (per figure-pixel moved)
SENS_TR  = 0.5;             % mm / pixel
SENS_ROT = 0.3;             % deg / pixel
SENS_SC  = 0.001;           % per pixel

%% ---- Orientation flip flags ----
% Identical to CT slicer — derived from MRI affine column vectors.
colX_ = Tmr(1:3,1);  colY_ = Tmr(1:3,2);  colZ_ = Tmr(1:3,3);
ax_xFlip  = colX_(1) > 0;
ax_yFlip  = colY_(2) < 0;
cor_xFlip = colX_(1) > 0;
cor_yTop  = colZ_(3) < 0;
sag_xFlip = colY_(2) < 0;

%% ---- Colour scheme ----
bg  = [0.10 0.10 0.10];
bg2 = [0.18 0.18 0.18];
bg3 = [0.22 0.22 0.22];
fg  = [0.92 0.92 0.92];

%% ---- Figure ----
ctCoregNii = '';    % default: cancelled
fig = figure('Name','CT-MR Registration', ...
    'NumberTitle','off', 'Color',bg, ...
    'Position',[30 30 1500 870], ...
    'WindowKeyPressFcn',   @cbKey, ...
    'WindowScrollWheelFcn',@cbScroll, ...
    'WindowButtonDownFcn', @cbDown, ...
    'WindowButtonMotionFcn',@cbMotion, ...
    'WindowButtonUpFcn',   @cbUp, ...
    'CloseRequestFcn',     @cbClose);

%% ---- Toolbar row 1: mode buttons + hint ----
tbH  = 0.048;
tb1Y = 0.950;
tb2Y = tb1Y - tbH - 0.006;

uicontrol(fig,'Style','text','String','Mode:', ...
    'Units','normalized','Position',[0.01 tb1Y 0.04 tbH], ...
    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',10);
btnTr = uicontrol(fig,'Style','togglebutton','String','Translate', ...
    'Units','normalized','Position',[0.055 tb1Y 0.082 tbH], ...
    'BackgroundColor',[0.25 0.55 0.25],'ForegroundColor','w', ...
    'FontSize',10,'Value',1,'Callback',@(~,~)setMode_('translate'));
btnRot = uicontrol(fig,'Style','togglebutton','String','Rotate', ...
    'Units','normalized','Position',[0.142 tb1Y 0.072 tbH], ...
    'BackgroundColor',bg2,'ForegroundColor',fg, ...
    'FontSize',10,'Value',0,'Callback',@(~,~)setMode_('rotate'));
btnSc = uicontrol(fig,'Style','togglebutton','String','Scale', ...
    'Units','normalized','Position',[0.219 tb1Y 0.065 tbH], ...
    'BackgroundColor',bg2,'ForegroundColor',fg, ...
    'FontSize',10,'Value',0,'Callback',@(~,~)setMode_('scale'));
uicontrol(fig,'Style','text', ...
    'String','Drag to transform CT  |  Click to move crosshair  |  Scroll: slice  |  Shift: fine', ...
    'Units','normalized','Position',[0.295 tb1Y 0.44 tbH], ...
    'BackgroundColor',bg,'ForegroundColor',[0.5 0.5 0.5],'FontSize',9, ...
    'HorizontalAlignment','left');

%% ---- Toolbar row 2: W/L + opacity ----
uicontrol(fig,'Style','text','String','MR W:', ...
    'Units','normalized','Position',[0.01 tb2Y 0.040 tbH], ...
    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
hMrW = uicontrol(fig,'Style','edit','String',num2str(mrWW), ...
    'Units','normalized','Position',[0.053 tb2Y 0.055 tbH], ...
    'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',9,'Callback',@cbWL);
uicontrol(fig,'Style','text','String','L:', ...
    'Units','normalized','Position',[0.110 tb2Y 0.020 tbH], ...
    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
hMrL = uicontrol(fig,'Style','edit','String',num2str(mrWL), ...
    'Units','normalized','Position',[0.133 tb2Y 0.055 tbH], ...
    'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',9,'Callback',@cbWL);

uicontrol(fig,'Style','text','String','CT W:', ...
    'Units','normalized','Position',[0.198 tb2Y 0.040 tbH], ...
    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
hCtW = uicontrol(fig,'Style','edit','String',num2str(ctWW), ...
    'Units','normalized','Position',[0.241 tb2Y 0.055 tbH], ...
    'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',9,'Callback',@cbWL);
uicontrol(fig,'Style','text','String','L:', ...
    'Units','normalized','Position',[0.298 tb2Y 0.020 tbH], ...
    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
hCtL = uicontrol(fig,'Style','edit','String',num2str(ctWL), ...
    'Units','normalized','Position',[0.321 tb2Y 0.055 tbH], ...
    'BackgroundColor',bg3,'ForegroundColor',fg,'FontSize',9,'Callback',@cbWL);

uicontrol(fig,'Style','text','String','CT opacity:', ...
    'Units','normalized','Position',[0.387 tb2Y 0.065 tbH], ...
    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',9);
hAlphaSlider = uicontrol(fig,'Style','slider','Min',0,'Max',1,'Value',ctAlpha, ...
    'Units','normalized','Position',[0.455 tb2Y+0.005 0.09 tbH*0.8], ...
    'Callback',@cbAlpha);
addlistener(hAlphaSlider,'Value','PostSet',@(~,~)cbAlphaLive());

%% ---- Slice panels ----
% Square panels (pixel-equal) so isotropic data fills without black bars.
% DataAspectRatio corrects for non-isotropic MRI voxels.
figPos = get(fig,'Position');
panW = 0.235;
panH = panW * figPos(3) / figPos(4);
panY = (tb2Y - panH) / 2 + 0.005;

axAx  = axes('Parent',fig,'Position',[0.005 panY panW panH], ...
    'Color','k','XColor','none','YColor','none'); hold(axAx,'on');
axCor = axes('Parent',fig,'Position',[0.248 panY panW panH], ...
    'Color','k','XColor','none','YColor','none'); hold(axCor,'on');
axSag = axes('Parent',fig,'Position',[0.491 panY panW panH], ...
    'Color','k','XColor','none','YColor','none'); hold(axSag,'on');

hImAx  = []; hImCor = []; hImSag = [];
hXHax  = []; hXHcor = []; hXHsag = [];
hTitleAx  = title(axAx,  '','Color',fg,'FontSize',9);
hTitleCor = title(axCor, '','Color',fg,'FontSize',9);
hTitleSag = title(axSag, '','Color',fg,'FontSize',9);

%% ---- Status bar ----
hStatus = uicontrol(fig,'Style','text', ...
    'Units','normalized','Position',[0.005 0.003 0.73 0.038], ...
    'BackgroundColor',bg,'ForegroundColor',[0.5 0.8 1.0], ...
    'FontSize',9,'HorizontalAlignment','left','String','');

%% ---- Right panel ----
rp = 0.742;  rw = 0.253;

uicontrol(fig,'Style','text','String','CT-MR Registration', ...
    'Units','normalized','Position',[rp 0.900 rw 0.048], ...
    'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',12,'FontWeight','bold');

labels_   = {'Tx (mm):','Ty (mm):','Tz (mm):','Rx (°):','Ry (°):','Rz (°):','Scale:'};
hReadout  = gobjects(1,7);
for ii = 1:7
    yy = 0.840 - (ii-1)*0.052;
    uicontrol(fig,'Style','text','String',labels_{ii}, ...
        'Units','normalized','Position',[rp yy 0.095 0.038], ...
        'BackgroundColor',bg,'ForegroundColor',fg,'FontSize',10, ...
        'HorizontalAlignment','right');
    hReadout(ii) = uicontrol(fig,'Style','text','String','0.0', ...
        'Units','normalized','Position',[rp+0.100 yy 0.080 0.038], ...
        'BackgroundColor',bg2,'ForegroundColor',[0.3 0.9 0.3],'FontSize',10, ...
        'HorizontalAlignment','center');
end

uicontrol(fig,'Style','pushbutton','String','Align Centers', ...
    'Units','normalized','Position',[rp 0.265 rw*0.95 0.048], ...
    'BackgroundColor',[0.3 0.45 0.6],'ForegroundColor','w', ...
    'FontSize',10,'Callback',@cbAlignCenters, ...
    'TooltipString','Snap CT centre to MRI centre (coarse alignment)');

uicontrol(fig,'Style','pushbutton','String','Reset Transform', ...
    'Units','normalized','Position',[rp 0.208 rw*0.95 0.048], ...
    'BackgroundColor',[0.55 0.25 0.25],'ForegroundColor','w', ...
    'FontSize',10,'Callback',@cbReset);

uicontrol(fig,'Style','pushbutton','String','Save', ...
    'Units','normalized','Position',[rp 0.100 rw*0.44 0.065], ...
    'BackgroundColor',[0.2 0.5 0.2],'ForegroundColor','w', ...
    'FontSize',12,'FontWeight','bold','Callback',@cbSave);
uicontrol(fig,'Style','pushbutton','String','Cancel', ...
    'Units','normalized','Position',[rp+rw*0.51 0.100 rw*0.44 0.065], ...
    'BackgroundColor',[0.5 0.25 0.25],'ForegroundColor','w', ...
    'FontSize',12,'Callback',@cbClose);

%% ---- Initial draw ----
refreshAll();
uiwait(fig);

%% ====================================================================
%% NESTED FUNCTIONS
%% ====================================================================

    % ------------------------------------------------------------------
    % Transform helpers
    % ------------------------------------------------------------------
    function T = buildRegTransform()
        % Row-vector convention: [ct_world 1] * T = [mr_world 1]
        % Applies: translate CT to origin → scale → rotate → translate back → user translation.
        ctr = ctCtrWorld(1:3);
        Tc  = eye(4);  Tc(4,1:3)  = -ctr;   % translate to origin
        Tci = eye(4);  Tci(4,1:3) =  ctr;   % translate back
        Ttr = eye(4);  Ttr(4,1:3) = [tx ty tz];
        S   = diag([sc sc sc 1]);

        cx = cosd(rx); sx_ = sind(rx);
        cy = cosd(ry); sy_ = sind(ry);
        cz = cosd(rz); sz_ = sind(rz);

        % Row-vector rotation matrices (post-multiply)
        Rx = [1  0    0   0;
              0  cx   sx_ 0;
              0 -sx_  cx  0;
              0  0    0   1];
        Ry = [cy  0  -sy_ 0;
              0   1   0   0;
              sy_ 0   cy  0;
              0   0   0   1];
        Rz = [cz  sz_ 0  0;
             -sz_ cz  0  0;
              0   0   1  0;
              0   0   0  1];

        T = Tc * S * Rx * Ry * Rz * Tci * Ttr;
    end

    function M = buildSamplingMatrix()
        % [mr_vox_1based, 1] * M = [ct_vox_1based, 1]
        % Chain: MR vox → MR world → inv(T_reg) → CT world → CT vox
        T_reg = buildRegTransform();
        M = Tmr / T_reg / Tct;   % = Tmr * inv(T_reg) * inv(Tct)
    end

    function ctSlice = sampleCTslice(dim, sliceIdx)
        % Sample CT onto an MRI slice plane via current transform.
        % dim: 1=axial(k=sliceIdx), 2=coronal(j=sliceIdx), 3=sagittal(i=sliceIdx)
        M = buildSamplingMatrix();
        switch dim
            case 1  % axial: i=1..nx, j=1..ny, k=fixed
                [jj, ii] = meshgrid(1:nyMr, 1:nxMr);
                kk = repmat(sliceIdx, nxMr, nyMr);
            case 2  % coronal: i=1..nx, k=1..nz, j=fixed
                [kk, ii] = meshgrid(1:nzMr, 1:nxMr);
                jj = repmat(sliceIdx, nxMr, nzMr);
            case 3  % sagittal: j=1..ny, k=1..nz, i=fixed
                [kk, jj] = meshgrid(1:nzMr, 1:nyMr);
                ii = repmat(sliceIdx, nyMr, nzMr);
        end
        sz_out = size(ii);
        N = numel(ii);
        mr_h = [ii(:), jj(:), kk(:), ones(N,1)];
        ct_h = mr_h * M;
        % interp3(V, Xq, Yq, Zq): X=col=dim2, Y=row=dim1, Z=slice=dim3
        ctSlice = reshape(interp3(ctVol, ct_h(:,2), ct_h(:,1), ct_h(:,3), ...
            'linear', NaN), sz_out);
    end

    function rgb = makeOverlay(mrSlice, ctSlice)
        mr_lo = mrWL - mrWW/2;  mr_hi = mrWL + mrWW/2;
        ct_lo = ctWL - ctWW/2;  ct_hi = ctWL + ctWW/2;

        mr_n = min(max((mrSlice - mr_lo) / (mr_hi - mr_lo), 0), 1);
        mr_rgb = repmat(mr_n, 1, 1, 3);

        if isempty(ctSlice) || all(isnan(ctSlice(:)))
            rgb = mr_rgb;
            return;
        end

        ct_n = min(max((ctSlice - ct_lo) / (ct_hi - ct_lo), 0), 1);
        hmap = hot(256);
        idx  = max(1, min(256, floor(ct_n * 255) + 1));
        ct_r = reshape(hmap(idx(:),1), size(ct_n));
        ct_g = reshape(hmap(idx(:),2), size(ct_n));
        ct_b = reshape(hmap(idx(:),3), size(ct_n));
        ct_rgb = cat(3, ct_r, ct_g, ct_b);

        % Alpha = ctAlpha * sqrt(ct_n): air/soft-tissue (ct_n≈0) stays transparent,
        % bone/metal (ct_n≈1) gets full ctAlpha. sqrt boosts mid-range bone brightness.
        a = ctAlpha * sqrt(ct_n) .* double(~isnan(ctSlice));
        a3 = repmat(a, 1, 1, 3);
        rgb = min(max((1-a3).*mr_rgb + a3.*ct_rgb, 0), 1);
    end

    % ------------------------------------------------------------------
    % Display
    % ------------------------------------------------------------------
    function refreshAll()
        if ~ishandle(fig), return; end

        % ---- Axial ----
        sl_mr = mrVol(:,:,curVox(3))';
        ct_raw = sampleCTslice(1, curVox(3))';
        if ax_yFlip, sl_mr = flipud(sl_mr); ct_raw = flipud(ct_raw); end
        rgb_ax = makeOverlay(sl_mr, ct_raw);
        if isempty(hImAx) || ~ishandle(hImAx)
            hImAx = image(axAx, rgb_ax);
            axis(axAx,'normal');
            set(axAx,'XLim',[0.5 nxMr+0.5],'YLim',[0.5 nyMr+0.5],'YDir','normal');
            set(axAx,'DataAspectRatio',[pdMr(1) pdMr(2) 1]);
            if ax_xFlip, set(axAx,'XDir','reverse'); end
            addOrientLabels(axAx,'R','L','A','P');
        else
            set(hImAx,'CData',rgb_ax);
        end
        set(hTitleAx,'String',sprintf('AXIAL  z=%d/%d',curVox(3),nzMr));

        % ---- Coronal ----
        sl_mr = squeeze(mrVol(:,curVox(2),:))';
        ct_raw = sampleCTslice(2, curVox(2))';
        if cor_yTop, sl_mr = flipud(sl_mr); ct_raw = flipud(ct_raw); end
        rgb_cor = makeOverlay(sl_mr, ct_raw);
        if isempty(hImCor) || ~ishandle(hImCor)
            hImCor = image(axCor, rgb_cor);
            axis(axCor,'normal');
            set(axCor,'XLim',[0.5 nxMr+0.5],'YLim',[0.5 nzMr+0.5],'YDir','normal');
            set(axCor,'DataAspectRatio',[pdMr(1) pdMr(3) 1]);
            if cor_xFlip, set(axCor,'XDir','reverse'); end
            addOrientLabels(axCor,'R','L','S','I');
        else
            set(hImCor,'CData',rgb_cor);
        end
        set(hTitleCor,'String',sprintf('CORONAL  y=%d/%d',curVox(2),nyMr));

        % ---- Sagittal ----
        sl_mr = squeeze(mrVol(curVox(1),:,:))';
        ct_raw = sampleCTslice(3, curVox(1))';
        if cor_yTop, sl_mr = flipud(sl_mr); ct_raw = flipud(ct_raw); end
        rgb_sag = makeOverlay(sl_mr, ct_raw);
        if isempty(hImSag) || ~ishandle(hImSag)
            hImSag = image(axSag, rgb_sag);
            axis(axSag,'normal');
            set(axSag,'XLim',[0.5 nyMr+0.5],'YLim',[0.5 nzMr+0.5],'YDir','normal');
            set(axSag,'DataAspectRatio',[pdMr(2) pdMr(3) 1]);
            if sag_xFlip, set(axSag,'XDir','reverse'); end
            addOrientLabels(axSag,'P','A','S','I');
        else
            set(hImSag,'CData',rgb_sag);
        end
        set(hTitleSag,'String',sprintf('SAGITTAL  x=%d/%d',curVox(1),nxMr));

        refreshCrosshairs();
        updateReadouts();
    end

    function refreshCrosshairs()
        if ~ishandle(fig), return; end
        xhColor = [0.3 1.0 0.3];
        lw = 0.8;

        iy_ax  = ternary(ax_yFlip,  nyMr - curVox(2) + 1, curVox(2));
        iz_cor = ternary(cor_yTop, nzMr - curVox(3) + 1, curVox(3));
        iz_sag = iz_cor;

        if ~isempty(hXHax)  && all(ishandle(hXHax)),  delete(hXHax);  end
        if ~isempty(hXHcor) && all(ishandle(hXHcor)), delete(hXHcor); end
        if ~isempty(hXHsag) && all(ishandle(hXHsag)), delete(hXHsag); end

        hXHax(1) = plot(axAx,[curVox(1) curVox(1)],[0.5 nyMr+0.5], ...
            'Color',xhColor,'LineWidth',lw,'HitTest','off');
        hXHax(2) = plot(axAx,[0.5 nxMr+0.5],[iy_ax iy_ax], ...
            'Color',xhColor,'LineWidth',lw,'HitTest','off');

        hXHcor(1) = plot(axCor,[curVox(1) curVox(1)],[0.5 nzMr+0.5], ...
            'Color',xhColor,'LineWidth',lw,'HitTest','off');
        hXHcor(2) = plot(axCor,[0.5 nxMr+0.5],[iz_cor iz_cor], ...
            'Color',xhColor,'LineWidth',lw,'HitTest','off');

        hXHsag(1) = plot(axSag,[curVox(2) curVox(2)],[0.5 nzMr+0.5], ...
            'Color',xhColor,'LineWidth',lw,'HitTest','off');
        hXHsag(2) = plot(axSag,[0.5 nyMr+0.5],[iz_sag iz_sag], ...
            'Color',xhColor,'LineWidth',lw,'HitTest','off');
    end

    function updateReadouts()
        vals  = {tx, ty, tz, rx, ry, rz, sc};
        fmts  = {'%.1f','%.1f','%.1f','%.1f','%.1f','%.1f','%.4f'};
        for ii = 1:7
            set(hReadout(ii),'String',sprintf(fmts{ii}, vals{ii}));
        end
        set(hStatus,'String', sprintf( ...
            '  T=[%.1f, %.1f, %.1f] mm   R=[%.1f, %.1f, %.1f]°   scale=%.4f', ...
            tx,ty,tz, rx,ry,rz, sc));
    end

    % ------------------------------------------------------------------
    % Mouse callbacks
    % ------------------------------------------------------------------
    function cbDown(~,~)
        ax = axUnderCursor();
        if isempty(ax), return; end
        fp = get(fig,'CurrentPoint');
        isDragging   = false;
        dragFigStart = fp;
        dragAxSel    = ax;
        dragState0   = [tx ty tz rx ry rz sc];
    end

    function cbMotion(~,~)
        if isempty(dragFigStart), return; end
        fp    = get(fig,'CurrentPoint');
        delta = fp - dragFigStart;        % [dx dy] in figure pixels
        if ~isDragging && norm(delta) < 3, return; end
        isDragging = true;

        fine  = any(strcmp(get(fig,'CurrentModifier'),'shift'));
        k = ternary(fine, 0.1, 1.0);

        dx = delta(1);  dy = delta(2);    % figure y increases upward → positive = drag up

        % Display-right always maps to -x world in axial/coronal (either because
        % colX_<0 makes dim1 go in -x, or XDir='reverse' flips it). Always negate.
        % Sagittal x (dim2/y) has the opposite convention — no negation needed.
        xSignAx  = -1;
        xSignSag =  1;

        ax = dragAxSel;
        switch mode_
            case 'translate'
                s = k * SENS_TR;
                if     ax == axAx,  tx = dragState0(1)+dx*s*xSignAx; ty = dragState0(2)+dy*s;
                elseif ax == axCor, tx = dragState0(1)+dx*s*xSignAx; tz = dragState0(3)+dy*s;
                elseif ax == axSag, ty = dragState0(2)+dx*s*xSignSag; tz = dragState0(3)+dy*s;
                end
            case 'rotate'
                s = k * SENS_ROT;
                if     ax == axAx,  rz = dragState0(6)+dx*s*(-xSignAx);
                elseif ax == axCor, ry = dragState0(5)+dx*s*xSignAx;
                elseif ax == axSag, rx = dragState0(4)+dx*s*(-xSignSag);
                end
            case 'scale'
                s = k * SENS_SC;
                sc = max(0.1, dragState0(7) + dx*s);
        end
        refreshAll();
    end

    function cbUp(~,~)
        if ~isDragging && ~isempty(dragFigStart) && ~isempty(dragAxSel)
            % Short click → move crosshair
            ax = dragAxSel;
            cp = get(ax,'CurrentPoint');
            xd = round(cp(1,1));  yd = round(cp(1,2));
            if ax == axAx
                xd = max(1,min(nxMr,xd));
                yd = max(1,min(nyMr,yd));
                curVox(1) = xd;
                curVox(2) = ternary(ax_yFlip, nyMr-yd+1, yd);
            elseif ax == axCor
                xd = max(1,min(nxMr,xd));
                yd = max(1,min(nzMr,yd));
                curVox(1) = xd;
                curVox(3) = ternary(cor_yTop, nzMr-yd+1, yd);
            elseif ax == axSag
                xd = max(1,min(nyMr,xd));
                yd = max(1,min(nzMr,yd));
                curVox(2) = xd;
                curVox(3) = ternary(cor_yTop, nzMr-yd+1, yd);
            end
            refreshAll();
        end
        isDragging   = false;
        dragFigStart = [];
        dragAxSel    = [];
        dragState0   = [];
    end

    function cbScroll(~,evt)
        ax = axUnderCursor();
        d  = evt.VerticalScrollCount;
        if     ax == axAx,  curVox(3) = max(1,min(nzMr, curVox(3)-d));
        elseif ax == axCor, curVox(2) = max(1,min(nyMr, curVox(2)-d));
        elseif ax == axSag, curVox(1) = max(1,min(nxMr, curVox(1)-d));
        end
        refreshAll();
    end

    function cbKey(~,evt)
        switch evt.Key
            case 'uparrow',    curVox(3) = min(nzMr, curVox(3)+1); refreshAll();
            case 'downarrow',  curVox(3) = max(1,    curVox(3)-1); refreshAll();
        end
    end

    % ------------------------------------------------------------------
    % Toolbar callbacks
    % ------------------------------------------------------------------
    function setMode_(m)
        mode_ = m;
        set(btnTr,  'Value', strcmp(m,'translate'), ...
            'BackgroundColor', ternary(strcmp(m,'translate'),[0.25 0.55 0.25],bg2));
        set(btnRot, 'Value', strcmp(m,'rotate'), ...
            'BackgroundColor', ternary(strcmp(m,'rotate'),[0.25 0.55 0.25],bg2));
        set(btnSc,  'Value', strcmp(m,'scale'), ...
            'BackgroundColor', ternary(strcmp(m,'scale'),[0.25 0.55 0.25],bg2));
    end

    function cbWL(~,~)
        mrWW = str2double(get(hMrW,'String'));
        mrWL = str2double(get(hMrL,'String'));
        ctWW = str2double(get(hCtW,'String'));
        ctWL = str2double(get(hCtL,'String'));
        refreshAll();
    end

    function cbAlpha(~,~)
        ctAlpha = get(hAlphaSlider,'Value');
        refreshAll();
    end

    function cbAlphaLive()
        ctAlpha = get(hAlphaSlider,'Value');
        refreshAll();
    end

    function cbReset(~,~)
        tx=0; ty=0; tz=0; rx=0; ry=0; rz=0; sc=1;
        refreshAll();
    end

    function cbAlignCenters(~,~)
        % Snap CT centre to MRI centre (coarse alignment starting point).
        mrCtrWorld = [round(nxMr/2), round(nyMr/2), round(nzMr/2), 1] * Tmr;
        tx = mrCtrWorld(1) - ctCtrWorld(1);
        ty = mrCtrWorld(2) - ctCtrWorld(2);
        tz = mrCtrWorld(3) - ctCtrWorld(3);
        refreshAll();
    end

    function cbSave(~,~)
        % Resample CT onto MRI grid using current transform, write NIfTI.
        set(hStatus,'String','  Resampling CT onto MRI grid... (this may take a few seconds)');
        drawnow;

        M  = buildSamplingMatrix();
        [ii, jj, kk] = ndgrid(1:nxMr, 1:nyMr, 1:nzMr);
        N  = numel(ii);
        mr_h  = [ii(:), jj(:), kk(:), ones(N,1)];
        ct_h  = mr_h * M;
        ctRs  = reshape(interp3(ctVol, ct_h(:,2), ct_h(:,1), ct_h(:,3), ...
            'linear', 0), nxMr, nyMr, nzMr);

        outPath  = fullfile(outputDir, 'ct_manual_coreg.nii');
        outInfo  = mrInfo;   % same grid as MRI
        niftiwrite(single(ctRs), outPath, outInfo, 'Compressed', false);

        ctCoregNii = outPath;
        fprintf('[ctMrRegistrationGUI] Saved: %s\n', outPath);
        uiresume(fig);
        delete(fig);
    end

    function cbClose(~,~)
        ctCoregNii = '';
        if ishandle(fig)
            uiresume(fig);
            delete(fig);
        end
    end

    % ------------------------------------------------------------------
    % Utility helpers
    % ------------------------------------------------------------------
    function ax = axUnderCursor()
        ax = [];
        fp  = get(fig,'CurrentPoint');   % figure pixels [x y]
        fsz = get(fig,'Position');       % [left bot w h] px
        np  = fp ./ fsz(3:4);           % normalized [x y]
        for aa = [axAx axCor axSag]
            apos = get(aa,'Position');   % [left bot w h] normalized
            if np(1)>=apos(1) && np(1)<=apos(1)+apos(3) && ...
               np(2)>=apos(2) && np(2)<=apos(2)+apos(4)
                ax = aa; return;
            end
        end
    end

    function addOrientLabels(ax, xl, xr, yt, yb)
        xl_ = get(ax,'XLim');  yl_ = get(ax,'YLim');
        mx = mean(xl_);  my = mean(yl_);
        text(ax, xl_(1)+0.01*(xl_(2)-xl_(1)), my, xl, ...
            'Color',[1 0.8 0],'FontSize',8,'HorizontalAlignment','left', ...
            'VerticalAlignment','middle','HitTest','off');
        text(ax, xl_(2)-0.01*(xl_(2)-xl_(1)), my, xr, ...
            'Color',[1 0.8 0],'FontSize',8,'HorizontalAlignment','right', ...
            'VerticalAlignment','middle','HitTest','off');
        text(ax, mx, yl_(2)-0.01*(yl_(2)-yl_(1)), yt, ...
            'Color',[1 0.8 0],'FontSize',8,'HorizontalAlignment','center', ...
            'VerticalAlignment','top','HitTest','off');
        text(ax, mx, yl_(1)+0.01*(yl_(2)-yl_(1)), yb, ...
            'Color',[1 0.8 0],'FontSize',8,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom','HitTest','off');
    end

end % ctMrRegistrationGUI

% ---- File-scope utility (not nested — can be called without capture) ----
function v = ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end
