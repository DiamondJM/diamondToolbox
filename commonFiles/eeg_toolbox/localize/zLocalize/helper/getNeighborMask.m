function neighborMask = getNeighborMask(elementInfo, myCoords, saveFile)
    % pairs/edges between close electrodes are found:
    %  + shear/bend/norm neighbors are considered pairs and found based on grid
    %    geometry.
    %  + inter-hardware edges are formed between electrodes within 15 mm
    %
    % Input:
    %   elementInfo - table with each row a hardware piece with columns colums:
    %   	tagName - unique value for each hardware piece
    %       bigDim  - largest dimension (>1 for grids, 1 for strips)
    %       smallDim- other dimension
    %       isCut   - string array like '[1 2]' where eval(isCut) gives indices
    %                 of missing electrodes on hardware row
    %   myCoords    - table with x, y, z, chanName columns
    %   saveFile (optional) - where to save edges variable
    %
    % Output: 
    %   - neighborMask is a somewhat more compacted form of the usable data.
    %       x=squareform(configVec) gives the usable, n by n data, such that x(i,j)
    %       is 1 if there is an edge between electrodes i and j and 0 otherwise.
    %       Without calling squareform, you can also use this as a mask for the
    %       vector output from the pdist function.

%  Copyright (C) 2017 Mike Trotta

    INTER_HARDWARE_THRESH = 15; % mm 

    % element based connectivity:
    % note the edges include duplicates
    
    [myCoords.tagName, myCoords.hwNdx] = util_split_stringnum(myCoords.chanName);
    nchan = height(myCoords);
    myCoords.chanNdx = [1 : nchan]';
    tagList = unique(myCoords.tagName);
    edges = [[] []];


    for iTag = 1:length(tagList)
        tag = tagList{iTag};


        % inter-hardware connectivity
        mask = strcmp(myCoords.tagName, tag);
        this_rows = myCoords(mask, :);
        other_rows= myCoords(~mask, :);
        this_xyz = [this_rows.x, this_rows.y, this_rows.z];
        other_xyz = [other_rows.x, other_rows.y, other_rows.z];

        dvec = pdist2(this_xyz, other_xyz);                                 
        [this_ndx, other_ndx] = find(dvec <= INTER_HARDWARE_THRESH);
        inter_edges = [this_rows.chanNdx(this_ndx), other_rows.chanNdx(other_ndx)];
        % ^ above index variables are index into chanNdx (all channels 1:n)

        % intra-hardware connectivity
        norm = findNeighbors([], [], tag, 1, 'elementInfo',elementInfo);
        bend = findNeighbors([], [], tag, 2, 'elementInfo',elementInfo);
        shear= findNeighbors([], [], tag, 1, 'elementInfo',elementInfo, 'useDiagonal',true);
        intra = [norm; bend; shear];
        len = length(intra);
        % ^ above variables are indices into the hardware
        
        temp_ndx = indexOf(intra, this_rows.hwNdx);
        % ^ this variable gives index into the hardware
        
        if any(temp_ndx(:) == 0)
            missing = unique(intra(find(temp_ndx == 0)));
            fprintf('The following numbers are missing: ');
            for i=1:numel(missing), fprintf('%d ',missing(i)); end
            fprintf('\n');
            error('The neighbors found based on element_info configuration were not found in CT_coordinates (%s). Check that the CT names match element_info',...
                tag);
        end
        
        intra_edges = this_rows.chanNdx(temp_ndx);
        % ^ above index variables are index into chanNdx (all channels 1:n)
        
        if size(intra_edges,2) == 1 && size(intra_edges, 1) == 2
            intra_edges = intra_edges';
        end
        
        edges = [edges; inter_edges; intra_edges];

    end

    % connect edges with upper triangular adjacency matrix vector 
    amat = zeros(nchan, nchan);
    if ~isempty(edges)
        for iEdge = 1:length(edges(:,1))
            amat(edges(iEdge,1), edges(iEdge,2)) = 1;
            amat(edges(iEdge,2), edges(iEdge,1)) = 1;
        end
    end
    neighborMask = squareform(amat,'tovector');

    % save edges to file
    if exist('saveFile','var') && ~isempty(saveFile)
        save(saveFile, 'edges');
    end

end

function [tagName, tagNdx] = util_split_stringnum(chanName)
    % splits chanName into tagName and electrode index in hardware

    chanName = cellstr(chanName);
    
    getDigitMask= @(s) isstrprop(s, 'digit');
    getDigits   = @(s,digitMask) str2double(s(digitMask));
    getName     = @(s,digitMask) s(1:find(digitMask,1)-1);
    
    digitMasks = getDigitMask(chanName);
    tagNdx     = cellfun(getDigits, chanName, digitMasks);
    tagName    = cellfun(getName, chanName, digitMasks, 'uniformOutput',0);
    
end
