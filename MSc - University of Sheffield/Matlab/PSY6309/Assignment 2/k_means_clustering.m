function [data_points_in_clusters_iterations, final_centers_coordinates_iterations] = k_means_clustering(data_points, k, manual_start_points)
    % Use random initial clusters centers if not given
    if (~exist('manual_start_points', 'var'))
        initial_centers_coordinates = zeros(k, 2);
        for i=1:numel(initial_centers_coordinates)
                initial_centers_coordinates(i) = rand(1);
        end
    else
        initial_centers_coordinates = manual_start_points;
    end
    
    % Set up parameters
    flag = 1;
    iterations = 0;
    last_centers_coordinates = initial_centers_coordinates;
    data_points_in_clusters_iterations = cell(1);
    final_centers_coordinates_iterations = cell(1);
    % Keep running until cluster centers stop changing
    while flag==1
        iterations = iterations + 1;
        [~, cluster_number] = assign_to_clusters(data_points, last_centers_coordinates);
        [data_points_in_clusters, updated_centers_coordinates] = re_calculate_centers_coordinates(data_points ,cluster_number, k);
        data_points_in_clusters_iterations{iterations} = data_points_in_clusters;
        final_centers_coordinates_iterations{iterations} = updated_centers_coordinates;
        if updated_centers_coordinates ~= last_centers_coordinates
            last_centers_coordinates = updated_centers_coordinates;
        else
            flag=0;
        end
    end
    
    %Visualize each iteration
    figure_rows = fix((iterations - 1)/3) + 1;
    figure()
    for iteration=1:iterations
        hold on;
        subplot(figure_rows, 3, iteration)
        visualize_k_means_clustering(data_points_in_clusters_iterations{iteration}, final_centers_coordinates_iterations{iteration}, iteration)
    end
end

function [distances_matrix, cluster_num] = assign_to_clusters(data_points, centers_coordinates)
    distances_matrix = zeros(size(centers_coordinates, 1), size(data_points, 1));
    cluster_num = zeros(size(data_points, 1), 1);
    for center_i=1:size(centers_coordinates, 1)
        for row_i=1:size(data_points, 1)
            distances_matrix(center_i, row_i) = pdist2(data_points(row_i, :), centers_coordinates(center_i, :), 'euclidean');
        end
    end
    for col_i=1:size(data_points, 1)
        [~, cluster_num(col_i)] = min(distances_matrix(:,col_i));
    end
end

function [data_points_in_clusters, new_centers_coordinates] = re_calculate_centers_coordinates(data_points, cluster_number, k)
    new_centers_coordinates = zeros(k, 2);
    data_points_in_clusters = cell(k, 1);
    for cluster_num=1:k
        cluster_points_indices = find(cluster_number==cluster_num);
        cluster_points = data_points(cluster_points_indices, :);
        new_centers_coordinates(cluster_num, :) = mean(cluster_points, 1);
        data_points_in_clusters{cluster_num} = cluster_points;
    end
end

function visualize_k_means_clustering(data_points_in_clusters, cluster_centers_coordinates, iteration_num)
    axis equal;
    %Plot Data Points in each cluster
    for cluster_num=1:numel(data_points_in_clusters)
        data_points_cluster_num = cell2mat(data_points_in_clusters(cluster_num));
        hold on;
        plot(data_points_cluster_num(:, 1), data_points_cluster_num(:, 2), 'o');
    end
    %Plot Cluster Centers
    hold on;
    plot(cluster_centers_coordinates(:, 1), cluster_centers_coordinates(:, 2), 'x');
    %Plot Seperating Lines
    hold on;
    try
        voronoi(cluster_centers_coordinates(:, 1), cluster_centers_coordinates(:, 2))
    catch 
    end
        title([num2str(iteration_num), ' Iteration(s)'])
    xlabel('PC 1')
    ylabel('PC 2')
end