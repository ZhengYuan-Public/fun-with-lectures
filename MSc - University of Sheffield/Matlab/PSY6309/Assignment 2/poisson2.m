function [spike_times_2 ,spike_train_2] = poisson2(z,T,dt)
    time_stamps = (0:dt:T);
    spike_train_2 = z*dt > rand(numel(time_stamps),1);
    spike_times_2 = find(spike_train_2==1)*dt;
end