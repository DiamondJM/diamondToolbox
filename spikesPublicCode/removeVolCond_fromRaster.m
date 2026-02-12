function sourceLoc = removeVolCond_fromRaster(sourceLoc)

rasters = sourceLoc.rasters;
waveforms = sourceLoc.waveforms;
listLocs = cell(1,size(rasters,2));

if sum(rasters(:)) == 0; return; end


volCondInds = find(sum(rasters,2) >= 2);
if isempty(volCondInds); return; end
fprintf('%.2f%% of single spike samples contained duplicate spikes.\n',full(length(volCondInds) / sum(any(rasters,2)) * 100));
% markForDeletion = cell(1,size(rasters,2));

for ii = 1:length(listLocs)
    listLocs{ii} = find(rasters(:,ii));
end


for ii = 1:length(listLocs)
    isBad = ismember(listLocs{ii},volCondInds);
    rasters(listLocs{ii}(isBad),ii) = false;
    waveforms{ii}(isBad,:) = [];
end

sourceLoc.rasters = rasters;
sourceLoc.waveforms = waveforms;



end