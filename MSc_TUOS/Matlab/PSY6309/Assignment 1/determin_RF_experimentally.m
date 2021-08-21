% Part II - Modelling receptive fields in visual cortex
% 2.3 Determine the neuron's receptive field experimentally
clear
clc
%{
Q1: The result is exactly the Gabor Filter generated from gabor_filter.m
Q2: This happens because these random images can be regarded as shapes in 
 all directions. The filter is a Gabor Filter, it's responses to
 images are already determined so when average the results we will still
 get the Gabor Filter.
%}
% Generate random images.
num_image = 10000;
image_row = 51;
image_col = 51;
random_images = rand(num_image, image_row * image_col);

% Prepare filter
sigma_x = 1;
sigma_y = 2;
k = 2;
[X, Y]=meshgrid(-5:.2:5,-5:.2:5);
filter = gabor_filter(X, Y, sigma_x, sigma_y, k);
filter = reshape(filter, 1, []);

% Calculate firing rate.
firing_rate = random_images .* filter;

% Multiply each stimulus with the respective firing rate then average.
temp = random_images .* firing_rate;
result = zeros(image_row, image_col);
for i=1:image_row * image_col
    result(i) = sum(temp(:,i)) / num_image;
end

% Plot the result
imagesc(result)