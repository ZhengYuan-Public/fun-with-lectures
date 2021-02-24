%% 1.5 Varying parameters of the LIF model
%
%% 1.5.1 Visualize
%
clear
synaptic_input_strength = [200, 400, 800, 1600];
tau_ss = [0.001, 0.003, 0.006, 0.012];

T = 2;
dt = 0.001;
thresh = 0.055;

N = numel(synaptic_input_strength);
input_spike_trains = {N};
for i=1:N
    z = synaptic_input_strength(i);
    input_spike_trains{i} = poisson2(z, T, dt);
end

psc_s = {N};
vm_s = {N};
output_spike_trains = {N};

for i=1:length(tau_ss)
    tau_s = tau_ss(i);
    for index=1:N
        [time_stamps, psc, vm, output_spike_train] = lif_modified(...
            input_spike_trains{index}, T, dt, tau_s);
        input_spike_train = input_spike_trains{index}; 
        psc_s{index} = psc;
        vm_s{index} = vm;
        output_spike_trains{index} = output_spike_train;

        visual_input_spike_train = input_spike_train;
        visual_psc = psc;
        visual_vm = vm;
        visual_output_spike_train = output_spike_train;

        lif_visualize(visual_input_spike_train, visual_psc, visual_vm, ...
            visual_output_spike_train, time_stamps, T, thresh);
        set(gcf, 'name', ...
            strcat('tas_s=', num2str(tau_s),...
            ', Input firing rate=', num2str(synaptic_input_strength(index)))...
            )
    end
end
%% 1.5.2 Describe and explain your results
%
% When the input spike train strength is small (firing rate = 200Hz), 
% the output spike is none, and as it becomes stronger (from 400Hz 
% to 800Hz and 1600Hz), more and more output spikes are generated. 
%
% With a fixed input strength, as the synaptic time constant changes
% from small (0.001) to large (0.012), more and more output spikes 
% are generated.
%
% As the input strength increases, the postsynaptic membrane potential
% is more and more likely to accumlate to a level greater than the 
% threshold, thus resunlting more output spikes. 
%
% Since the synaptic time constant characterizes the synaptic 
% conductance, the larger the value is inputs conducted to the 
% postsynaptic membrane are greater, resulting more output spikes.
