%% 1.3 The inhomogeneous Poisson process
% 1.3.1 Test your function with different arrays z and visualize the results
clear
T = 2;
dt = 0.001;
time_stamps = (0:dt:T);
zs_list = {...
    100*randn(numel(time_stamps)), ...
    abs(100*sin(time_stamps)), ...
    abs(100*cos(time_stamps)), ...
    100*tan(time_stamps), ...
    100*1./time_stamps, ...
    100*log(time_stamps+1), ...
    };
for zs_index=1:size(zs_list, 2)
    zs = zs_list{zs_index};
    poisson3(zs, T, dt);
end
%%
% 1.3.2 Switch from a baseline firing rate (100 Hz) to a high firing rate (800 Hz) in the middle of a 2 second-long spike train.
clear
T = 2;
dt = 0.001;
zs = zeros(T/dt+1, 1);
for index=1:numel(zs)
    if index < numel(zs)/2
        zs(index) = 100;
    else
        zs(index) = 800;
    end
end
[spike_times_3 ,spike_train_3] = poisson3(zs, T, dt);