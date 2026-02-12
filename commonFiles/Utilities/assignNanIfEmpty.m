function val = assignNanIfEmpty(val)

if isempty(val); val = nan; end

end