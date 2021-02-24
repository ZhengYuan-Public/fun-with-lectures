function [spike_times ,binned_spike_train] = poisson1(z,T,dt)
    time_stamps = (0:dt:T);
    binned_spike_train = zeros(numel(time_stamps),1);
    time_rvs = rand(10*numel(time_stamps), 1);
    intervals = zeros(10*numel(time_stamps),1);
    for time = 1:numel(time_rvs)
        intervals(time) = -log(time_rvs(time))/z;
    end
    % Only use intervals within time T
    valid_num_intervals = find(cumsum(intervals) < T, 1, 'last' );
    spike_times = zeros(valid_num_intervals, 1);
    for time_index = 1:valid_num_intervals
        spike_times(time_index) = sum(intervals(1:time_index));
    end
    for spike_index=1:numel(spike_times)
        num_bin = ceil(spike_times(spike_index) / dt);
        binned_spike_train(num_bin) = 1;
    end
end