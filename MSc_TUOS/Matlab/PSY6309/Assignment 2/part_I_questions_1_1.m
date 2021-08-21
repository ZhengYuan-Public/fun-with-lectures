%% 1.1 Generating Poisson spike trains
% 1.1.1 Test functions with different parameters
clear
T = 2; % second 
firing_rate_values = [10, 50, 100, 200, 400, 800]; % Hz
dt_values = [0.05, 0.01, 0.005, 0.001]; % s

for firing_rate_index=1:length(firing_rate_values)
    fr = firing_rate_values(firing_rate_index);
    for dt_values_index=1:length(dt_values)
        dt = dt_values(dt_values_index);
        [spike_times ,binned_spike_train] = poisson1(fr, T, dt);
        raster_plot(spike_times, binned_spike_train, fr, dt);
        set(gca, 'xlim', [0, T]);
    end
end
%%
% 1.1.2 How to choose dt
% 
% Based on the obeservation of different values of dt, 0.001s (1ms)
% should be an appropriate dt.
%