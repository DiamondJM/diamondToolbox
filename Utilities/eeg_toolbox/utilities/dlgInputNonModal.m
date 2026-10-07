function answers = dlgInputNonModal(message, title, labels, defaults)
% dlgInputNonModal  Non-modal text-input dialog in the dlgNonModal theme.
%
% Replacement for inputdlg that matches dlgNonModal (dark theme, no
% WindowStyle='modal', so no macOS click-through).
%
% Usage:
%   answers = dlgInputNonModal(message, title, labels, defaults)
%
% message  - char string or cell array of strings shown above the fields
% title    - window title string, also shown as the dialog heading
% labels   - cell array of field labels
% defaults - cell array of default field strings (same length as labels)
%
% Returns a cell array of the entered strings, or {} if cancelled/closed.

if ischar(message), message = {message}; end
if ischar(labels),  labels  = {labels};  end
if nargin < 4 || isempty(defaults), defaults = repmat({''}, size(labels)); end
nF = numel(labels);

BG    = [0.13 0.13 0.13];
FG    = [0.92 0.92 0.92];
BTN   = [0.25 0.25 0.25];
GREEN = [0.18 0.42 0.18];
FIELD = [0.20 0.20 0.20];

W = 560;  PAD = 12;  HDR_H = 30;  BTN_H = 34;  ROW_H = 28;  LBL_W = 80;
FS = 11;

ss  = get(0, 'ScreenSize');
fig = figure('Name', title, 'MenuBar', 'none', 'ToolBar', 'none', ...
    'NumberTitle', 'off', 'Color', BG, 'Resize', 'off', 'Visible', 'off', ...
    'Position', [round((ss(3)-W)/2) round(ss(4)/2) W 200], ...
    'CloseRequestFcn', @cbClose);

% message height from wrapped line count
probe = uicontrol(fig, 'Style', 'text', 'String', 'Ag', 'FontSize', FS, ...
    'Units', 'pixels', 'Position', [0 0 W 100]);
lineH = get(probe, 'Extent');  lineH = lineH(4);
set(probe, 'Position', [0 0 W-PAD*2 100]);
wrapped = textwrap(probe, message);
delete(probe);
MSG_H = ceil(numel(wrapped) * lineH * 1.25) + 8;

fieldsH = nF * (ROW_H + PAD/2);
totalH  = PAD + BTN_H + PAD + fieldsH + PAD/2 + MSG_H + PAD + HDR_H + PAD;
set(fig, 'Position', [round((ss(3)-W)/2) round((ss(4)-totalH)/2) W totalH]);

btnY = PAD;
y    = btnY + BTN_H + PAD;            % bottom of lowest field row
hEdit = gobjects(1, nF);
for i = nF:-1:1
    uicontrol(fig, 'Style', 'text', 'String', labels{i}, 'FontSize', FS, ...
        'ForegroundColor', FG, 'BackgroundColor', BG, ...
        'HorizontalAlignment', 'left', ...
        'Position', [PAD y LBL_W ROW_H-6]);
    hEdit(i) = uicontrol(fig, 'Style', 'edit', 'String', defaults{i}, ...
        'FontSize', FS, 'ForegroundColor', FG, 'BackgroundColor', FIELD, ...
        'HorizontalAlignment', 'left', ...
        'Position', [PAD+LBL_W y W-PAD*2-LBL_W ROW_H]);
    y = y + ROW_H + PAD/2;
end

msgY = y + PAD/2;
uicontrol(fig, 'Style', 'text', 'String', message, 'FontSize', FS, ...
    'ForegroundColor', FG, 'BackgroundColor', BG, ...
    'HorizontalAlignment', 'left', ...
    'Position', [PAD msgY W-PAD*2 MSG_H]);
uicontrol(fig, 'Style', 'text', 'String', title, 'FontSize', 13, ...
    'FontWeight', 'bold', 'ForegroundColor', FG, 'BackgroundColor', BG, ...
    'HorizontalAlignment', 'left', ...
    'Position', [PAD msgY+MSG_H+PAD W-PAD*2 HDR_H]);

bW = 110;
uicontrol(fig, 'Style', 'pushbutton', 'String', 'OK', 'FontSize', FS, ...
    'FontWeight', 'bold', 'ForegroundColor', FG, 'BackgroundColor', GREEN, ...
    'Position', [W-PAD-bW*2-PAD btnY bW BTN_H], 'Callback', @cbOK);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', 'FontSize', FS, ...
    'ForegroundColor', FG, 'BackgroundColor', BTN, ...
    'Position', [W-PAD-bW btnY bW BTN_H], 'Callback', @cbClose);

set(fig, 'Visible', 'on');
uicontrol(hEdit(1));

setappdata(0, 'dlgInputNonModal_result', {});
uiwait(fig);
answers = getappdata(0, 'dlgInputNonModal_result');
if isappdata(0, 'dlgInputNonModal_result')
    rmappdata(0, 'dlgInputNonModal_result');
end

    function cbOK(~,~)
        setappdata(0, 'dlgInputNonModal_result', get(hEdit, 'String'));
        delete(fig);
    end

    function cbClose(~,~)
        setappdata(0, 'dlgInputNonModal_result', {});
        delete(fig);
    end
end
