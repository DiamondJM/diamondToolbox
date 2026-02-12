function warnOnce(msg)

if isequal(lastwarn,msg); return; end
warning(msg); 

end