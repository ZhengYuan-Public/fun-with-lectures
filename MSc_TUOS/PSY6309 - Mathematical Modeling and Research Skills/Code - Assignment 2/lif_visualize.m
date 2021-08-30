function lif_visualize(input_spike_train, psc, vm, output_spike_train, time_stamps, T, thresh)
    figure()
    hold on;
    
    subplot(4, 1, 1);
    p1 = plot([input_spike_train, input_spike_train], [0, 0.5], 'r');
    xlim([0, T]);
    
    subplot(4, 1, 2);
    p2 = plot(time_stamps, psc, 'g');
    xlim([0, T]);
    
    subplot(4, 1, 3);
    p3 = plot(time_stamps, vm, 'b');
    xlim([0, T]);
    if ~isempty(output_spike_train)
        hold on;
        for i = 1:length(output_spike_train)
            plot([output_spike_train', output_spike_train'],[thresh, thresh + 0.1],'Color','k'); 
        end
    end
    
    if ~isempty(output_spike_train)
        subplot(4, 1, 4);
        p4 = plot([output_spike_train', output_spike_train'], [0, 0.5], 'm');
        xlim([0, T]);
    else
        subplot(4, 1, 4);
        p4 = plot(0, 0, 'm'); % empty plot
        xlim([0, T]);
    end
    
    legend([p1(1), p2(1), p3(1), p4(1)], ...
        'Input spike train',...
        'Postsynaptic current',...
        'Membrane potential',...
        'Output spike train');  
end