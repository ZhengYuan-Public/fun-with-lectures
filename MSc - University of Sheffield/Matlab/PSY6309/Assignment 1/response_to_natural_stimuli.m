% Part II - Modelling receptive fields in visual cortex
% 2.4 Modelling the responses to natural stimuli
clear
clc
%{
Q1: The result shows the edges of trees having vertical lines.
Q2: The result looks like this because the kernel used here respondes most
intensively to vertical shapes.
%}

% Load images, set kernel size and calculate steps for moving the kernel.
load('forest.mat')
kernel_row = 51;
kernel_col = 51;
[pix_row, pix_col] = size(im);
steps_row = pix_row - kernel_row + 1;
steps_col = pix_col - kernel_col + 1;

% Prepare kernel/filter.
sigma_x = 1;
sigma_y = 2;
k = 2;
[X, Y]=meshgrid(-5:.2:5,-5:.2:5);
filter = gabor_filter(X, Y, sigma_x, sigma_y, k);

result = zeros(steps_row, steps_col); % Preallocation.
for row=1:steps_row
    for col=1:steps_col
        result(row, col) = ...
            sum(sum(im(row:row+kernel_row-1, ...
            col:col+kernel_col-1) .* filter)) / (kernel_row*kernel_col) ;
    end
end

result = rectify(result);
imagesc(result);