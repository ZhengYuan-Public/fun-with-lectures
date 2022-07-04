%% Parameters
t_start = 0; % In second;
t_end = 1; % In second;
dt = 0.001; % Step size;
t = (t_start : dt : t_end); % Discretized time;
%% LIF Neuron with about 10Hz firing rate
V_rest1 = -0.07; % Resting potential (V) of neuron 1;
V_threshold1 = -0.04; % Firing threshold (V) of neuron 1;
R1=100E6; % Resistance (ohm) of neuron 1;
tau_m1 = 0.02; % Membrane time constant of neuron 1;
I_inj = 0.305E-9; % Injected current (A);

V_m1=zeros(1, numel(t));
V_m1(1) = V_rest1; % Boundary condition, start from V_rest1;
spike_timings_index1 = []; % Exact spike timings index of neuron 1;
for i=2:numel(t)
    V_m1(i) = V_m1(i-1)+dt/tau_m1*(V_rest1-V_m1(i-1)+R1*I_inj);
    if V_m1(i)>V_threshold1
        spike_timings_index1(end+1) = i; % Record the time as spike time;
        V_m1(i) = V_rest1-0.01; % Rset after spike, smaller than V_rest because of hyperpolarization;
    end
end
%% Conductance-based synapse
t_f = t(spike_timings_index1);
g_syn_bar = 40E-12; % In second, =40ps;
tau_syn = 5E-3; % In second, =5ms;
g_syn = zeros(1, numel(t));
for i=1:numel(t)
    spike_conuts = find(t(i) > t_f);
    temp = 0;
    for spike_num=1:numel(spike_conuts)
        temp = g_syn_bar*exp(-(t(i)-t_f(spike_num))/tau_syn)*heaviside(t(i)-t_f(spike_num)) + temp;
    end
        g_syn(i) = temp;
end
%% Post-synaptic neuron
V_rest2 = -0.070; % Resting potential (V) of neuron 2;
V_threshold2 = -0.05; % Firing threshold (V) of neuron 2;
R2 = 300E9; % Resistance (ohm) of neuron 2;
tau_m2 = 0.08; % Membrane time constant of neuron 2;

E_syn = -0.1;
I_t =g_syn.*(V_m1-E_syn); % Post-synaptic current (A);

V_m2=zeros(1, numel(t));
V_m2(1) = V_rest2; % Boundary condition, start from V_rest2
spike_timings_index2 = []; % Exact spike timings index for neuron 2
for i=2:numel(t)
    V_m2(i) = V_m2(i-1)+dt/tau_m2*(V_rest2-V_m2(i-1)+R2*I_t(i));
    if V_m2(i)>V_threshold2
        spike_timings_index2(end+1) = i; %#ok<*SAGROW> % Record the time as spike time
        V_m2(i) = V_rest2 - 0.005; % Rset after spike, smaller than V_rest because of hyperpolarization;
    end
end
%% Plot
% V_m1
figure(1)
hold on
plot(t, V_m1, 'b', 'linewidth', 1)
title('Pre-synaptic LIF Neuron')
x_s1 = t(spike_timings_index1)-dt;
y1 = V_threshold1;
y2 = y1+0.07;
plot([x_s1', x_s1'], [y1, y2], 'b', 'linewidth', 1)
ylabel('V(t)')
ylim([-0.09 0.05])
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'q2_A3_1.png')

% V_m2
figure(2)
hold on
plot(t, V_m2, 'b', 'linewidth', 1)
title('Post-synaptic LIF Neuron')
x_s2 = t(spike_timings_index2)-dt;
y1 = V_threshold2-0.02;
y2 = y1+0.07;
plot([x_s2', x_s2'], [y1, y2], 'b', 'linewidth', 1)
ylabel('V(t)')
ylim([-0.09 0.05])
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'q2_A3_2.png')

% g_syn
figure(3)
hold on
plot(t, g_syn, 'b', 'linewidth', 1)
title('g_{syn}')
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'q2_A3_3.png')

% reaster plot
figure(4)
hold on
y1 = 0;
y2 = y1 + 0.01;
y3 = y2 + 0.001;
y4 = y3 + 0.01;
pre = plot([x_s1(:), x_s1(:)], [y3, y4], 'b', 'linewidth', 1, 'DisplayName', 'Pre-synaptic Neuron');
post = plot([x_s2(:), x_s2(:)], [y1, y2], 'r', 'linewidth', 1, 'DisplayName', 'Post-synaptic Neuron');
title('Raster plot')
legend([pre(1), post(1)], 'Pre-synaptic Neuron', 'Post-synaptic Neuron');
set(gca,'ytick',[])
ylim([0 0.028])
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'q2_A3_4.png')