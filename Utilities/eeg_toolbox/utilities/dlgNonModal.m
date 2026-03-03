function result = dlgNonModal(message, title, varargin)
% dlgNonModal  Non-modal question/info dialog (avoids macOS click-through).
%
% Replacement for questdlg that does not set WindowStyle = 'modal', so
% clicking on the MATLAB desktop behind the dialog does not pass through
% to windows beneath MATLAB.
%
% Usage:
%   result = dlgNonModal(message, title, btn1)
%   result = dlgNonModal(message, title, btn1, btn2)
%   result = dlgNonModal(message, title, btn1, btn2, btn3)
%
% message  - char string or cell array of strings (one element per line)
% title    - window title string
% btn1..3  - button label strings; the leftmost button is the default
%
% Returns the label of the button clicked, or '' if the window is closed.

buttons = varargin;
nBtn    = numel(buttons);
assert(nBtn >= 1 && nBtn <= 3, 'dlgNonModal: provide 1–3 button labels.');

if ischar(message), message = {message}; end

% ---- theme ---------------------------------------------------------------
BG  = get(0, 'defaultUicontrolBackgroundColor');
FG  = [0.1 0.1 0.1];

% ---- layout --------------------------------------------------------------
W      = 420;
PAD    = 14;
BTN_H  = 28;
BTN_W  = 110;
LINE_H = 18;
nLines = numel(message);
MSG_H  = nLines * LINE_H + PAD;
totalH = PAD + MSG_H + PAD + BTN_H + PAD;

ss  = get(0, 'ScreenSize');
fig = figure( ...
    'Name',        title, ...
    'MenuBar',     'none', ...
    'ToolBar',     'none', ...
    'NumberTitle', 'off', ...
    'Color',       BG, ...
    'Resize',      'off', ...
    'Position',    [round((ss(3)-W)/2) round((ss(4)-totalH)/2) W totalH], ...
    'CloseRequestFcn', @cbClose);

% ---- message text --------------------------------------------------------
uicontrol(fig, 'Style', 'text', ...
    'String',              message, ...
    'ForegroundColor',     FG, ...
    'BackgroundColor',     BG, ...
    'FontSize',            11, ...
    'HorizontalAlignment', 'left', ...
    'Position',            [PAD PAD+BTN_H+PAD W-PAD*2 MSG_H]);

% ---- buttons -------------------------------------------------------------
% Distribute buttons right-to-left so btn1 is on the left (default/primary).
btnX = W - PAD - BTN_W;
for i = nBtn:-1:1
    bi = i;
    uicontrol(fig, 'Style', 'pushbutton', ...
        'String',   buttons{bi}, ...
        'FontSize', 11, ...
        'Position', [btnX PAD BTN_W BTN_H], ...
        'Callback', @(~,~) cbBtn(buttons{bi}));
    btnX = btnX - BTN_W - 6;
end

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
