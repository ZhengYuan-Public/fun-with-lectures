clear
close all;
%% Virus-antibodies
dt = 1e-3; % Step size
t_all = 120; % Total time of simulation (hour)
times = 0 : dt : t_all; % Discretized time

N = 10e5; % Threshold value of virus that patient dies
n_0 = 6; % Initial condition for virus
m_0 = 0; % Initial condition for anti-bodies

alpha  = 0.515; beta = 0.17; gamma = 0.979; delta = 0.3; % Parameters

n = zeros(1, numel(times)); % The number of virus
n(1) = n_0; % Initial condition, n(t=0) = n_0
m = zeros(1, numel(times)); % the number of anti-bodies
m(1) = m_0; % Initial condition, m(t=0) = m_0
for t=2:numel(times)
    n(t) = n(t-1)+dt * (alpha*n(t-1)-beta*m(t-1)*n(t-1)) * heaviside(N-n(t-1)); 
    % Forward Euler method
    m(t) = m(t-1)+dt * (gamma*n(t-1)-delta*m(t-1)*n(t-1)) * heaviside(N-n(t-1)); 
    % Forward Euler method
end
%% Plot the result
figure(1); hold on;
p1 = plot(times, n, 'r', 'linewidth', 2);
p2 = plot(times, m, 'g', 'linewidth', 2);
legend({'Virus', 'Anti-bodies'})
title({'Virus and Anti-bodies population'})
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'q2c2.png')