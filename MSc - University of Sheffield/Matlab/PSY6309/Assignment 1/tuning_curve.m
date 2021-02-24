% Part II - Modelling receptive fields in visual cortex
% 2.2 Plot the neuron's tuning curve
clear
clc
% Load images.
load('edges.mat')
images = reshape(ms, 72, []);

% Prepare filter
sigma_x = 1;
sigma_y = 2;
k = 2;
[X, Y]=meshgrid(-5:.2:5,-5:.2:5);
filter = gabor_filter(X, Y, sigma_x, sigma_y, k);
filter = reshape(filter, 1, []);

% Calculate firing rate
temp = images .* filter;
firing_rate = zeros(72, 1);
for i=1:72
    firing_rate(i) = sum(temp(i,:));
end

plot(theta, firing_rate);