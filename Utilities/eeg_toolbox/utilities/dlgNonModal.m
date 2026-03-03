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
% btn1..3  - button labels; the leftmost is the primary/default
%
% Returns the label of the button clicked, or '' if the window is closed.

buttons = varargin;
nBtn    = numel(buttons);
assert(nBtn >= 1 && nBtn <= 3, 'dlgNonModal: provide 1-3 button labels.');

if ischar(message), message = {message}; end

% Join lines for the edit box
msgStr = strjoin(message, char(10));

% ---- layout --------------------------------------------------------------
W      = 400;
PAD    = 12;
BTN_H  = 26;
BTN_W  = 100;
MSG_H  = max(80, min(240, numel(message) * 16 + 20));
totalH = PAD + MSG_H + PAD + BTN_H + PAD;

ss  = get(0, 'ScreenSize');
fig = figure( ...
    'Name',        title, ...
    'MenuBar',     'none', ...
    'ToolBar',     'none', ...
    'NumberTitle', 'off', ...
    'Resize',      'off', ...
    'Position',    [round((ss(3)-W)/2) round((ss(4)-totalH)/2) W totalH], ...
    'CloseRequestFcn', @cbClose);

% ---- message (read-only edit box — handles wrapping and overflow) --------
uicontrol(fig, 'Style', 'edit', ...
    'String',              msgStr, ...
    'FontSize',            10, ...
    'HorizontalAlignment', 'left', ...
    'Max',                 10, ...
    'Min',                 0, ...
    'Enable',              'inactive', ...
    'Position',            [PAD PAD+BTN_H+PAD W-PAD*2 MSG_H]);

% ---- buttons (right-aligned, btn1 leftmost) ------------------------------
btnX = W - PAD - BTN_W;
for i = nBtn:-1:1
    bi = i;
    uicontrol(fig, 'Style', 'pushbutton', ...
        'String',   buttons{bi}, ...
        'FontSize', 10, ...
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
