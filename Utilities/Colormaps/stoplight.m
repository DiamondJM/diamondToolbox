function cmap = stoplight(n)

% if nargin == 0; n = 256; end 
% assert(iseven(n)); 
% 
% cmap = [interp1([0;1],[1 0 0; 1 1 1],linspace(0,1,n / 2));interp1([0;1],[1 1 1; 0 1 0],linspace(0,1,n / 2))];

% check the inputs
if nargin == 0
    n = 256;
elseif nargin > 1
    error('too many inputs')
elseif ~isnumeric(n) || n <=3
    error('input must be an integer greater than 3')
end

% set the x, y, and z valules and interpolate
x           = [1,ceil(n/2),n];
xnew        = 1:n;
y           = 1:3;
ynew        = 1:3;
[x,y]       = meshgrid(x,y);
[xnew,ynew] = meshgrid(xnew,ynew);
z(:,1)      = [1 0 0]';
z(:,2)      = [1 1 0]';
z(:,3)      = [0 1 0]';
cmap        = interp2(x,y,z,xnew,ynew)';


end