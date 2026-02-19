function jitter = simpleBar(x,y,varargin)

p = inputParser;
addParameter(p,'errorType','dots');
addParameter(p,'selectDot',0);
addParameter(p,'color','b'); 
parse(p,varargin{:})
errorType = p.Results.errorType;
selectDot = p.Results.selectDot; 
color = p.Results.color; 

hold on
meanY = mean(y,'omitnan');
bar(x,meanY,'faceColor',color); 

if isequal(errorType,'bar')
    err = calculateConfInterval(y);
    % err = std(y,'omitnan'); 
    err = errorbar(x,meanY,err); err.LineStyle = 'none'; err.Color = 'k'; 
elseif isequal(errorType,'dots')
    jitter = randn(size(y))/20;
    plot(x + jitter(1:length(y) ~= selectDot), y(1:length(y) ~= selectDot),'k.','MarkerSize',20);
    
    if selectDot; plot(x + jitter(selectDot),y(selectDot),'kx','MarkerSize',10); end
    
end

end
