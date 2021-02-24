%% 1.2 An alternative method for generating Poisson spike trains
% 1.2.1 Visualize poisson2.m
clear
T = 2; % second 
firing_rate_values = [10, 50, 100]; % Hz
dt_values = [0.01, 0.005, 0.001]; % s

for firing_rate_index=1:length(firing_rate_values)
    fr = firing_rate_values(firing_rate_index);
    for dt_values_index=1:length(dt_values)
        dt = dt_values(dt_values_index);
        [spike_times ,binned_spike_train] = poisson2(fr, T, dt);
        raster_plot(spike_times, binned_spike_train, fr, dt)
    end
end
%%
% 1.2.2 Compare the two methods for different z and dt and describe similarities, differences, advantages, and disadvantages of the two methods.
clear
T = 2; % second 
firing_rate_values = [50, 100, 200]; % Hz
dt_values = [0.05, 0.01, 0.001]; % s
for firing_rate_index=1:length(firing_rate_values)
    fr = firing_rate_values(firing_rate_index);
    for dt_values_index=1:length(dt_values)
        dt = dt_values(dt_values_index);
        [spike_times ,binned_spike_train] = poisson1(fr, T, dt);
        [spike_times_2, ~] = poisson2(fr,T,dt);
        raster_plot(spike_times, binned_spike_train, fr, dt, spike_times_2)
    end 
end

%%
% 1. Similarity and Differences of these two methods: 
% When dt is small, these two methods have similar perfomance; 
% while dt is larger, the inversion method works fine 
% but the other one has poor performance.
%
% 2. Advantages and Disadvantages of these two methods: 
% The inversion method is computationally expensive, but it can
% handle larger dt compared to the other one; while the other method is
% relatively simple but cannot handle large dt.
%
