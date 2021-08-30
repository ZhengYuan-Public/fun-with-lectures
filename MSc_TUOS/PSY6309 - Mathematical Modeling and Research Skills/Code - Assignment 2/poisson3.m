function [spike_times_3 ,spike_train_3] = poisson3(zs, T, dt)
    time_stamps = (0:dt:T);
    thr_s = (zs*dt)';
    rv_s = rand(numel(time_stamps), 1);
    spike_train_3 = zeros(numel(time_stamps), 1);
    for index = 1:numel(spike_train_3)
        if thr_s(index) > rv_s(index)
            spike_train_3(index) = 1;
        end
    end
    spike_times_3 = find(spike_train_3 == 1)*dt;
    simple_raster_plot(spike_times_3 ,spike_train_3, dt)
    set(gca, 'xlim', [0, T])
end

function simple_raster_plot(spike_times, binned_spike_train, dt)
    figure()
    hold on;
    width=1920;
    height=100;
    left=0;
    bottem=0;
    set(gcf,'position', [left,bottem,width,height])
    xlabel('Time(s)')
    set(gca,'ytick',[]);
    set(gca,'TickLength',[0 0])
    binned_spike_times = find(binned_spike_train == 1)*dt;
    
    p1 = plot([spike_times, spike_times], [0, 0.4], 'b', 'LineWidth', dt/3);
    p2 = plot([binned_spike_times, binned_spike_times], [0.5, 0.9], 'r', 'LineWidth', dt/3);
    legend([p1(1), p2(1)], 'Spike times', 'Binned spike train')
end