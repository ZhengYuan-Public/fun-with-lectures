%% 2.2 A more challenging spike sorting problem
%
clear;
load('waveforms_difficult.mat')

% Run PCA
cov_matrix = cov(waveforms);
[eigenvector_matrix, diagnal_matrix] = eig(cov_matrix);
[sorted_diagnal_matrix, sorted_diagnal_matrix_indices] = sort(diag(diagnal_matrix),'descend');
selected_dim_indices = find((cumsum(sorted_diagnal_matrix)/sum(sorted_diagnal_matrix))>0.84,1);
selected_eigenvector_matrix = eigenvector_matrix(:,sorted_diagnal_matrix_indices(1:selected_dim_indices));
projected_data = waveforms*selected_eigenvector_matrix; 

[idx,C,sumd,D] = kmeans(projected_data, 3);

%% 
% 2.2.1 Decide how many principal components to keep and justify your chioce
%  
% Each value of info_percent represent information that each principle
% component (PC) contains. Start from the 4th PC, each PC contains less than 5.5%
% of the total information, thus I decide to choose 3 PCs, the cummulated information
% sum(info_percent(1:3)) = 0.8417. When visualzing the projected_data, We
% can observe some clustering spatial structures.
%
%% 
% 2.2.2 How does the complexity of this data set differ from the one used above?
%
% In this data set, information are more evenly distributed in more dimensions, so 
% we need another principle components to represent the data. So compared to the last
% date set, this one is more complicated.
%

%% 
% 2.2.3 How many clusters can you identify?
% 
% 3 clusters can be identified
%

%% 
% 2.2.4 Visualize the result
%
figure()
scatter3(projected_data(:,1), projected_data(:,2), projected_data(:,3), 15, idx, 'filled');

%% 
% 2.2.5 Ways to estimate the number of clusters
% 
% There are different ways to choose the optimal number of clusters k, such
% as the elbow method and the silhouette. In MATLAB, we can use
% evalclusters() to evaluate the clustering result of each k then choose
% the optimal k.
%

%% 
% 2.2.6 Pitfalls and problems with the k-means method
%
% 1. Clustering data of varying sizes and density: k-means has trouble 
% clustering data where clusters are of varying sizes and density.
%
% 2. Dependent on initial values: for a low k, you can mitigate this 
% dependence by running k-means several times with different initial 
% values and picking the best result. As k increases, you need advanced 
% versions of k-means to pick better values of the initial centroids (called k-means seeding).
%
% 3. Clustering outliers: centroids can be dragged by outliers, 
% or outliers might get their own cluster instead of being ignored.
%