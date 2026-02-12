function [xunit,yunit] = circle(x,y,r,th,plotting)

hold on
if nargin <= 4; plotting = false; end
if nargin <= 3; th = [0 2 * pi]; end

th = linspace(th(1),th(2),100); 
% th = 0 : pi / 50 : 2 * pi; 
xunit = r * cos(th) + x;
yunit = r * sin(th) + y; 

% fill(xunit,yunit,[1 0 0]);

if ~plotting; return ;end
plot(xunit,yunit,'k'); 

end
