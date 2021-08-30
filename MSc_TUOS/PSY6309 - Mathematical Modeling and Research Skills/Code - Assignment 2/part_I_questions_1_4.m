%% 1.4 Poisson input to the leaky integrate-and-fire model
%%
% 1.4.1 Briefly explain how synaptic input is incorporated in this implementation
% In lif.m, synaptic inputs are given by postsynaptic current (PSC) events
% happening times "T_psc". Then these times are fit into the predefined, 
% evenly distributed time points "times", which can be accessed by its 
% index "index_pscs". Finally, iterate through all time points "times" 
% and calculate the membrane potentials v_m. If v_m is greater than a 
% predefined threshold value, an output spike is generated at the time.
%

%%
% 1.4.2Visualization
%
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
close all;

z = 1000;
tau_s = 0.003;
thresh = 0.055;
N = 3; % number of trials
[spike_times_1, ~] = poisson1(z, T, dt);
[spike_times_2, ~] = poisson1(z, T, dt);
input_spike_trains = {spike_times_1, spike_times_2, spike_times_3};

psc_s = {N};
vm_s = {N};
output_spike_trains = {N};


for index=1:N
    [time_stamps, psc, vm, output_spike_train] = lif_modified(input_spike_trains{index}, T, dt);
    input_spike_train = input_spike_trains{index}; 
    psc_s{index} = psc;
    vm_s{index} = vm;
    output_spike_trains{index} = output_spike_train;
    
    visual_input_spike_train = input_spike_train;
    visual_psc = psc;
    visual_vm = vm;
    visual_output_spike_train = output_spike_train;
    
    lif_visualize(visual_input_spike_train, visual_psc, visual_vm, visual_output_spike_train,...
        time_stamps, T, thresh); 
end