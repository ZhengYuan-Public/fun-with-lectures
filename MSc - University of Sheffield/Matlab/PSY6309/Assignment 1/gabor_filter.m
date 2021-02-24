% Part II - Modelling receptive fields in visual cortex
% 2.1 Generating the Gabor filter
%{
>>[X, Y]=meshgrid(-5:.2:5,-5:.2:5);
>>imagesc(gabor_filter(X, Y, 1, 2, 2));
%}
function z = gabor_filter(x, y, sigma_x, sigma_y, k)
    z = gauss2D(x, y, sigma_x, sigma_y) .* cosine2D(k, x);
end