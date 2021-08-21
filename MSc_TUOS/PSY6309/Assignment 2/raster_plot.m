function raster_plot(spike_times, binned_spike_train, z, dt, spike_times_2)
    figure()
    hold on;
    width=1920;
    height=100;
    left=0;
    bottem=0;
    set(gcf,'position', [left,bottem,width,height])
    xlabel('Time(s)')
    title(['Firing rate=', num2str(z), ', dt=', num2str(dt)])
    set(gca,'ytick',[]);
    set(gca,'TickLength',[0 0])
    binned_spike_times = find(binned_spike_train == 1)*dt;
    
    if (~exist('spike_times_2', 'var'))
        p1 = plot([spike_times, spike_times], [0, 0.4], 'b', 'LineWidth', dt/3);
        p2 = plot([binned_spike_times, binned_spike_times], [0.5, 0.9], 'r', 'LineWidth', dt/3);
        legend([p1(1), p2(1)], 'Spike times', 'Binned spike train')
    else
        p1 = plot([spike_times, spike_times], [0, 0.4], 'b', 'LineWidth', dt/3);
        p2 = plot([binned_spike_times, binned_spike_times], [0.5, 0.9], 'r', 'LineWidth', dt/3);
        p3 = plot([spike_times_2, spike_times_2], [1, 1.4], 'g', 'LineWidth', dt/3);
        legend([p1(1), p2(1), p3(1)], 'Spike times', 'Binned spike train', 'Spike times 2');
    end
end