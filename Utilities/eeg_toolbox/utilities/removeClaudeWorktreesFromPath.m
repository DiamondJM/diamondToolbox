function removeClaudeWorktreesFromPath()
% removeClaudeWorktreesFromPath  Drop .claude/ folders from the MATLAB path.
%
% addpath(genpath(<toolbox root>)) also adds .claude/worktrees/*, which hold
% full (possibly stale) copies of the toolbox.  Because '.claude' sorts
% before 'Utilities', those copies shadow the real Utilities functions.
% Called from the electrodeLocalizer and sourceLocalizer constructors.

parts = strsplit(path, pathsep);
bad   = parts(contains(parts, [filesep '.claude' filesep]) | ...
              endsWith(parts, [filesep '.claude']));
if ~isempty(bad)
    rmpath(bad{:});
    fprintf('[path] Removed %d .claude folder(s) from the MATLAB path (stale worktree copies).\n', numel(bad));
end
end
