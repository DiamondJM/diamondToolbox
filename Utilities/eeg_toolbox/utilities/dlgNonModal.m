function result = dlgNonModal(message, title, varargin)
% dlgNonModal  Non-modal question/info dialog (avoids macOS click-through).
%
% Replacement for questdlg that does not set WindowStyle = 'modal', so
% clicking on the MATLAB desktop behind the dialog does not pass through
% to windows beneath MATLAB.  Dark-themed to match the electrodeLocalizer
% setup dialog; the message area and buttons size themselves to their
% text, so nothing is clipped.
%
% Usage:
%   result = dlgNonModal(message, title, btn1)
%   result = dlgNonModal(message, title, btn1, btn2)
%   result = dlgNonModal(message, title, btn1, btn2, btn3)
%
% message  - char string or cell array of strings (one element per line;
%            long lines wrap)
% title    - window title string, also shown as the dialog heading
% btn1..3  - button labels; the leftmost is the primary/default
%
% Returns the label of the button clicked, or '' if the window is closed.

buttons = varargin;
nBtn    = numel(buttons);
assert(nBtn >= 1 && nBtn <= 3, 'dlgNonModal: provide 1-3 button labels.');

if ischar(message), message = {message}; end

% ---- theme (matches electrodeLocalizer.localizationSetupDialog) -----------
BG    = [0.13 0.13 0.13];
FG    = [0.92 0.92 0.92];
BTN   = [0.25 0.25 0.25];
GREEN = [0.18 0.42 0.18];

% ---- layout --------------------------------------------------------------
W      = 560;
PAD    = 12;
HDR_H  = 30;
BTN_H  = 34;
MSG_FS = 11;

ss  = get(0, 'ScreenSize');
fig = figure( ...
    'Name',        title, ...
    'MenuBar',     'none', ...
    'ToolBar',     'none', ...
    'NumberTitle', 'off', ...
    'Color',       BG, ...
    'Resize',      'off', ...
    'Visible',     'off', ...
    'Position',    [round((ss(3)-W)/2) round(ss(4)/2) W 200], ...
    'CloseRequestFcn', @cbClose);

% ---- message: wrap to the dialog width, then size the box to fit ----------
hMsg = uicontrol(fig, 'Style', 'text', ...
    'String',              message, ...
    'FontSize',            MSG_FS, ...
    'ForegroundColor',     FG, ...
    'BackgroundColor',     BG, ...
    'HorizontalAlignment', 'left', ...
    'Units',               'pixels', ...
    'Position',            [PAD 0 W-PAD*2 100]);
wrapped = textwrap(hMsg, message);
set(hMsg, 'String', wrapped);
ext   = get(hMsg, 'Extent');
MSG_H = min(ext(4) + 6, round(ss(4) * 0.6));

totalH = PAD + BTN_H + PAD*2 + MSG_H + PAD + HDR_H + PAD;
btnY   = PAD;
msgY   = btnY + BTN_H + PAD*2;
hdrY   = msgY + MSG_H + PAD;
set(hMsg, 'Position', [PAD msgY W-PAD*2 MSG_H]);
set(fig,  'Position', [round((ss(3)-W)/2) round((ss(4)-totalH)/2) W totalH]);

% ---- heading -------------------------------------------------------------
uicontrol(fig, 'Style', 'text', ...
    'String',              title, ...
    'FontSize',            13, ...
    'FontWeight',          'bold', ...
    'ForegroundColor',     FG, ...
    'BackgroundColor',     BG, ...
    'HorizontalAlignment', 'left', ...
    'Position',            [PAD hdrY W-PAD*2 HDR_H]);

% ---- buttons: sized to label, right-aligned, btn1 leftmost and green -----
hBtn = gobjects(1, nBtn);
bW   = zeros(1, nBtn);
for i = 1:nBtn
    isPrimary = (i == 1);
    hBtn(i) = uicontrol(fig, 'Style', 'pushbutton', ...
        'String',          buttons{i}, ...
        'FontSize',        11, ...
        'FontWeight',      ternary(isPrimary, 'bold', 'normal'), ...
        'ForegroundColor', FG, ...
        'BackgroundColor', ternary(isPrimary, GREEN, BTN), ...
        'Callback',        @(~,~) cbBtn(buttons{i}));
    e = get(hBtn(i), 'Extent');
    bW(i) = max(110, e(3) + 32);
end
btnX = W - PAD - sum(bW) - PAD*(nBtn-1);
for i = 1:nBtn
    set(hBtn(i), 'Position', [btnX btnY bW(i) BTN_H]);
    btnX = btnX + bW(i) + PAD;
end
uicontrol(hBtn(1));   % keyboard focus on the default action

set(fig, 'Visible', 'on');

% ---- block ---------------------------------------------------------------
setappdata(0, 'dlgNonModal_result', '');
uiwait(fig);
result = getappdata(0, 'dlgNonModal_result');
if isappdata(0, 'dlgNonModal_result')
    rmappdata(0, 'dlgNonModal_result');
end

    function cbBtn(label)
        setappdata(0, 'dlgNonModal_result', label);
        delete(fig);
    end

    function cbClose(~,~)
        setappdata(0, 'dlgNonModal_result', '');
        delete(fig);
    end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
