% Part I - Functions
% 1.3 Creating a two-dimensional Gaussian bump
%{
>>[X, Y]=meshgrid(-5:.2:5,-5:.2:5);
>>imagesc(gauss2D(X, Y, 1, 2));
These codes generate an image with a set of ovals sharing the
same center. It's yellow at the center and gradually changes to blue.
 
By increasing/decrecing X or Y, the image will stretch/shrink in the
corespounding axis;
By increasing/decrecing sigma_x or sigma_y, the ovals will stretch/shrink
in the corespounding axis;

By changeing the input range, the image changes size accordingly, but the 
coordinates of the center does not change.

%}
function z = gauss2D(x, y, sigma_x, sigma_y)
    z = 1/(2*sigma_x*sigma_y*pi)*exp(-x.^2/(2*sigma_x.^2)-y.^2/(2*sigma_y.^2));
end