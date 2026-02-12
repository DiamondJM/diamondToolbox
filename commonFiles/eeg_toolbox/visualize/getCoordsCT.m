function [xyz, chanNames, table] = getCoordsCT(subj, rootEEGdir)
% Get a subject's CT-registered electrode locations
%
% USAGE: [xyz, chanNames, t] = getCoordsCT(subj, rootEEGdir);
%
% DESCRIPTION
%   Returns n xyz coordinates and their n channel names
%
% OUTPUT
%   xyz - n xyz coords
%   chanName - n chanNames corresponding to xyz
%   table - table with chanName, x, y, and z columns

    tal = fullfile(rootEEGdir, subj, 'tal');
    filename = fullfile(tal, 'intermediates/locs_1_registered/coords.mr_pre.ict2bmrAF.csv');
    if ~exist(filename,'file')
        error('Coordinate file not found here: %s\n', filename);
    end
    
    t = readtable(filename);
    xyz = [t.x t.y t.z];
    chanNames = t.chanName;
    table = t(:, 1:4);
    
end