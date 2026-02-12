function colormap2D = create2dColormap

fontSize = 20;
%==================================================================================================
% Show lab color space for values of L = 50.
numColors = 256;
ramp = linspace(-100,100, numColors);
cform = makecform('lab2srgb');
a = repmat(ramp, [numColors 1]);           % -a on left
b = repmat(flipud(ramp'), [1 numColors]);  % -b on bottom
L = 50 * ones(numColors, numColors);  % A single L value.
Lab = cat(3, L, a, b); % A 2D image.
colormap2D = applycform(Lab, cform);

return; 


% Display it.
subplot(2, 1, 1);
imshow(colormap2D);
axis on;
caption = sprintf('2D Colormap');
title(caption, 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
drawnow;
% Now that RGB image is made up (to be a 2-D colormap), make scatterplot.
numPoints = 1000;
amplitude = 3000;
x = amplitude * rand(1, numPoints);
y = amplitude * rand(1, numPoints);
subplot(2, 1, 2);
grid on;
thisColor = zeros(numPoints, 3);
for k = 1 : numPoints
	col = ceil(x(k) * numColors / amplitude);
	row = ceil(y(k) * numColors / amplitude);
	thisColor(k, :) = [colormap2D(row, col, 1), colormap2D(row, col, 2), colormap2D(row, col, 3)];
	fprintf('(x,y) = (%6.1f, %6.1f), row = %3d, col = %3d, thisColor = (%.4f, %.4f, %.4f)\n', x(k), y(k), row, col, thisColor(k, :));
% 	plot(x(k), y(k), '.', 'Color', thisColor, 'MarkerSize', 30);
% 	hold on;
end
scatter(x, y, 30 * ones(1, numPoints), thisColor, 'filled');
axis('square');
grid on;
xlim([0, amplitude]);
