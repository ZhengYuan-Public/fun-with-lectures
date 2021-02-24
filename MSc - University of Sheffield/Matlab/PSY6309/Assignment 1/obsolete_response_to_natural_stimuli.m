clear
clc
% 2.4 Modelling the responses to natural stimuli
load('forest.mat')
chunck_row = 51;
chunck_col = 51;
[pix_row, pix_col] = size(im);
num_row = fix(pix_row / chunck_row);
num_col = fix(pix_col / chunck_col);
new_im = im(1:51*num_row, 1:51*num_col); % Delete extra pixels.

% Split IM into small images.
images = mat2cell(new_im, ...
    chunck_row*ones(num_row, 1), chunck_col*ones(num_col, 1));

% Prepare filter
sigma_x = 1;
sigma_y = 2;
k = 2;
[X, Y]=meshgrid(-5:.2:5,-5:.2:5);
filter = gabor_filter(X, Y, sigma_x, sigma_y, k);

% Create firing_rates matrix correspounding to new_im.
firing_rates = zeros(num_row*chunck_row, num_col*chunck_col);
firing_rates = mat2cell(firing_rates, ...
    chunck_row*ones(num_row, 1), chunck_col*ones(num_col, 1));

for index=1:numel(images)
    firing_rates{index} = images{index} .* filter;
end

result = rectify(cell2mat(firing_rates));
% imagesc(new_im);
imagesc(result);