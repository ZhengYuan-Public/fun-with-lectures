function [times, I_syn, V, t_spikes] = lif_modified(spike_times, T, dt, tau_s)
    %%% set up parameters (SI units)

    % membrane constants
    
    tau_m = 0.030;
    Rm = 3e7; % membrane resistance

    % Time variables
    delta_t = dt; % length of single simulation timestep (in seconds)
    duration = T; % duration of the simulation
    No_steps = round(duration ./ delta_t);
    times = linspace(0, duration, No_steps + 1); % time points for the update

    % membrane potential
    V_0 = 0.05; % The initial membrane potential 
    V = zeros(1, No_steps + 1);
    V(1) = V_0;

    thresh = 0.055; % voltage threshold for spike generation
    if ~exist('tau_s', 'var') % synaptic time constant
        tau_s = 0.003; 
    end
    %%% set up PSC events

    % Times of PSCs

    T_psc = spike_times;
    Q = 40e-13; % charge deliverd by a single PSC event
    I_0 = Q ./ tau_s; % resulting normalisation constant for the exponential description of the PSC
    index_pscs = round(T_psc ./ delta_t); 
    I_syn = zeros(1, No_steps + 1); % initialise synaptic current

    t_spikes = [];
    for i = 1:No_steps
        if ismember(i, index_pscs)
            I_syn(i) = I_syn(i) + I_0;
        end

        % update the synaptic current
        I_syn(i+1) = I_syn(i) - (1 ./ tau_s) .* I_syn(i) .* delta_t;

        dV = (1 ./ tau_m) .* (-V(i) + Rm .* I_syn(i)) .* delta_t; 
        V(i+1) = V(i) + dV;
        if V(i+1) > thresh
            V(i+1) = 0;
            t_spikes = [t_spikes (i - 1) * delta_t]; 
        end
    end
end