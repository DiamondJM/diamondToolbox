function [xyz, index] = getLocation(subj, rootEEGdir, chanName, varargin)
% GETLOCATION returns the location of a patient's given channel ("lead")
%
% INPUTS:
%   subj
%   rootEEGdir
%   chanName - example FG1
%
% OPTIONAL INPUT:
%   surf - Optional. If you have one, pass it to avoid loading it.
%   leads - Optional. If you have them (phys), pass to avoid loading them
%   fast - Optional. Skips cross validation AND index output for speed.
%
% OUTPUTS
%   xyz - [x,y,z] of point in subject's space (will be on the dura, approximately)
%   index - index or indices in subject's std 141 mesh. This is a integer
%       [1,198812], and corresponds to whichever hemisphere the lead is in
%
% DESCRIPTION
%
% Gets the xyz-coordinate of a lead. For Utahs, the nearest mesh index is
% returned. For larger electrodes, a set of nearest mesh indices are
% returned. 
%
% NOTE: If you need this function to be faster and you don't need the
% index, you can modify it to skip loading the surface/finding the indices.
% You can also skip cross validation. Right now, this is built to look at
% monopolar.mat as well as xyz-coordinates to cross validate our
% localization results between sources which should be identical.
%
% REVISION HISTORY:
%   09/16 MST - Created
%   10/16 MST - Add fast param

ip = inputParser;
ip.addParameter('surf', []);
ip.addParameter('fast', false);
ip.addParameter('leads', []);
ip.parse(varargin{:});
fast = ip.Results.fast;
surf = ip.Results.surf;
leads = ip.Results.leads;

if isempty(leads)
    leads = getLeads(subj, rootEEGdir);
end
tal = fullfile(rootEEGdir, subj, 'tal');
fnameMonopolar = fullfile(tal, 'monopolar.mat');
fnameLeads = fullfile(tal, 'leads.csv');
fnameLeadsBackup = fullfile(tal, 'intermediates/locs_4_compiled/compiled.csv');

assert(~isempty(find(strcmp(chanName, leads), 1)), 'chanName not found in leads');

coords = []; %#ok<*NASGU>
xyz = [];
index = [];
xyz_monopolar = [];
ndx_monopolar = [];

% try monopolar (as a double-check)
if ~fast && exist(fnameMonopolar, 'file')
    monopolar = getfield(load(fnameMonopolar), 'monopolar');
    rowMask = strcmpi({monopolar.chanName}, chanName);
    nMatches = length(find(rowMask));
    if nMatches == 1
        xyz_monopolar = monopolar(rowMask).native_xyz;
        ndx_monopolar = monopolar(rowMask).std_nodes;
        hemi = monopolar(rowMask).loc_hemi;
    elseif nMatches > 1
        error('SEVERE ERROR: more than 1 row match in monopolar.mat');
    end
end
% try primary
if exist(fnameLeads, 'file')
    coords = readtableSafe(fnameLeads);
    if ~isempty(coords)
        row = coords(strcmpi(coords.chanName, chanName), :);
        if ~isempty(row)
            xyz = [row.x, row.y, row.z];
            if fast, return; end
        end
    end
end
% try secondary on failure
if isempty(xyz) && exist(fnameLeadsBackup, 'file')
    fprintf('Could not find %s in %s. \nLooking in %s now\n', chanName, fnameLeads, fnameLeadsBackup);
    coords = readtableSafe(fnameLeadsBackup);
    if ~isempty(coords)
        row = coords(strcmpi(coords.chanName, chanName), :);
        if ~isempty(row)
            xyz = [row.x, row.y, row.z];
        end
    end
end

% cross validate
if ~isempty(xyz) && ~isempty(xyz_monopolar)
    xyzFmt = '[%d, %d, %d]';
    msg = ['Cross validation of leads xyz and monopolar failed:\n'...
        xyzFmt, ' ~= \n'...
        xyzFmt];
    assert(isequal(xyz, xyz_monopolar(1,:)), msg, xyz, xyz_monopolar);
end

% in case we ONLY have monopolar...
if isempty(xyz)
    xyz = xyz_monopolar;
end

if isempty(xyz)
    error('Could not get location of %s', chanName);
end
    
% == We have the XYZ, now get the nearest node(s) by looking at surface ==
sumaDir = fullfile(rootEEGdir, subj, 'tal/raw/SUMA');
if ~exist(sumaDir, 'dir')
    sumaDir = fullfile(rootEEGdir, subj, 'tal/intermediates/images_2_fsSumaStd', subj, 'SUMA');
end
        
% lookup jacktable
jacktable = getJackTable(subj, rootEEGdir);
if isempty(jacktable), return; end

% get hemisphere and hardware type
row = jacktable(strcmpi(chanName,jacktable.chanName),:);
type = row.hardwareType;
hemi = char(row.whichHemi);

if ~exist('surf', 'var')
    surf = surf_load_main('whichHemi', hemi, 'dir', sumaDir);
end

if strcmpi(type, 'utah')
    [index, dist] = knnsearch(surf.vertices, xyz);
    if isempty(ndx_monopolar)
        fprintf('Vertex found %d mm away', dist);
    else
        % cross validate
        assert(isequal(ndx_monopolar, index), ...
            ['Index in monopolar does not match nearest xyz vertex:\n'...
            'ndx %d found %d mm away verses\n'...
            'monopolar ndxs: ' repmat('%d ', 1, length(ndx_monopolar))],...
            index, dist, ndx_monopolar);
            
    end
    
else
    index = ndx_monopolar;
end

end % end function getLocation