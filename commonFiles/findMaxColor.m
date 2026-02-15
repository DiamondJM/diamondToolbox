function maxColor = findMaxColor(bp, VPerROI, valPerROIC, varargin)
% [vertexData, vertexCnt] = plotRegionsData(bp, VPerROI, valPerROIC, varargin)
% plots data on multiple vertex sets
%
% Input:
%   VPerROI     - n x 1 cell array of surface index arrays; Vertices for each ROI
%   valPerROIC  - n x 1 array of data values to plot; Value for each ROI center
%
% Optional Key-value pairs:
%   surf        - Surface name to plot on
%   rh_begin    - If you're giving [lh; rh], begin rh at this index
%   clim        - [low high]
%   cmap        - colormap
%   blend       - [default=1] 1 gives full color, 0.1 gives 10 percent color
%   cumulative  - if 0 [default], take mean of overlapping vertex sets
%                 if 1, take sum of overlapping vertex sets
%                 if 2, take max
%                 if 3, median
% Output:
%   vertexData  - data plotted at each vertex from [1:max(vertexSets)]
%   vertexCnt   - how many data points are at each vertex
%
% Details:
%   Uses plotRegion to plot a new surface and uses colormap-indexed colors
%
%



ip = inputParser;
ip.addParameter('surf',[]);
ip.addParameter('rh_begin',[]);
ip.addParameter('clim',get(gca,'CLim'));
ip.addParameter('cmap',colormap);
ip.addParameter('blend',1);
ip.addParameter('cumulative',0);
ip.addParameter('haxis', []);
ip.parse(varargin{:});

surfName = ip.Results.surf;
rh_begin = ip.Results.rh_begin;
clim = ip.Results.clim;
cmap = ip.Results.cmap;
blend = ip.Results.blend;
cumulative = ip.Results.cumulative;
haxis = ip.Results.haxis;
if isempty(haxis)
    haxis = gca();
end
% Check for lh/rh split
if rh_begin > 1
    % myBp, VPerROI, valPerROI,'surf','rh'
    
    modargin = {'surf','lh','rh_begin',0};
    maxColorL = findMaxColor(bp, VPerROI(1:rh_begin-1), valPerROIC(1:rh_begin-1), modargin{:});
    modargin = {'surf','rh','rh_begin',0};
    maxColorR = findMaxColor(bp, VPerROI(rh_begin:end), valPerROIC(rh_begin:end), modargin{:});
%     vertexData = cat(1, vertexData_l, vertexData_r);
%     vertexCnt = cat(1, vertexCnt_l, vertexCnt_r);
    
    maxColor = max([maxColorL maxColorR]);
    return;
end

% main
n = length(VPerROI);
assert(length(valPerROIC) == n, 'data and VPerROI must be the same size');
%surf = bp.get_psurf(surfName);


vertexData = zeros(bp.stdNumVert, 1);
vertexCnt  = zeros(bp.stdNumVert, 1);
vertexDataC= cell(bp.stdNumVert, 1);


% aggregate data from different vertex sets onto each vertex
for i = 1:n
    vertexCnt(VPerROI{i}) = vertexCnt(VPerROI{i}) + 1;
    if cumulative == 2
        % Max
        vertexData(VPerROI{i}) = nanmax(vertexData(VPerROI{i}), valPerROIC(i));
    elseif cumulative == 3
        % Cells
        error('Not implemented: median (cumulative==3)')
        %vertexDataC(VPerROI{i}) = [vertexDataC(VPerROI{i}) repmat({vertexData(VPerROI{i})}, numel(VPerROI{i},1)];
    else
        % Sum (for mean, divide later)
        vertexData(VPerROI{i}) = vertexData(VPerROI{i}) + valPerROIC(i);
    end
end

% average at each vertex
findMask = vertexCnt > 0;
vertexData(~findMask) = nan;
if cumulative == 0 % average the sum
    vertexData(findMask) = vertexData(findMask) ./ vertexCnt(findMask);
elseif cumulative == 3
    %vertexData(findMask) = cellfun(@nanmedian,
    error('Not implemented')
end

maxColor = max(vertexData);

end
