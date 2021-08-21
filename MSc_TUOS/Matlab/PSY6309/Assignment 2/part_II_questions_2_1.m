%% 2.1 Implement the k-means algorithm for spike sorting
%

clear;
clc;
load('waveforms.mat')

% PCA
cov_matrix = cov(waveforms);
[eigenvector_matrix, diagnal_matrix] = eig(cov_matrix);
[sorted_diagnal_matrix, sorted_diagnal_matrix_indices] = sort(diag(diagnal_matrix),'descend');
selected_dim_indices = find((cumsum(sorted_diagnal_matrix)/sum(sorted_diagnal_matrix))>0.809,1);
selected_eigenvector_matrix = eigenvector_matrix(:,sorted_diagnal_matrix_indices(1:selected_dim_indices));
projected_data = waveforms*selected_eigenvector_matrix; 

%%
% 2.1.1 Explore how sensitive the method is to: 
%
%%
% 1. The initial cluster centres
%
for i=1:4
    k_means_clustering(projected_data, 3); % Random intial cluster centers
    set(gcf, 'name', 'k-means clustering. K=3')
end

%%
% 2. The number of clusters used
%
for cluster=3:6
    k_means_clustering(projected_data, cluster);
    set(gcf, 'name', strcat('Example of k-means clustering. K=', num2str(cluster)));
end
%%
% 2.1.2 How well does your method perform in sorting spikes?
% 
% After running the algorithms for about 100 times, for most of the time
% the program is able to divide the data into three clusters within 6
% iterations, but in some cases it failed to reach desired result.
% 
% Here is an example when the inital cluster centers are intialed to 
% [-0.5827, 0.3728; 1.6651, -0.2981; 2.0743, 0.2497].
%
manual_start_points = [-0.5827, 0.3728; 1.6651, -0.2981; 2.0743, 0.2497];
k_means_clustering(projected_data, 3, manual_start_points);
set(gcf, 'name', 'An example of anomaly k-means clustering. K=3')
