clear
close all
%% Parameters
tau = 1;
delta_t = 0.001;
t = (0: delta_t: 10*tau); % Discretized time;
%% Analytical Method
% Boundary condition: tau = 1; V(0)=1;
V_analytical= zeros(1, numel(t));
for i=1:numel(t)
    N = min(fix(t(i)/tau), 10);
    V_analytical(i) = ((1-exp(N+1))/(1-exp(1))/tau+0)*exp(-t(i)/tau);
end
%% Euler's Method
V_euler = zeros(1, numel(t));
V_euler(1) = 1; % Same boundary condition as the analytical method;
delta_func_index = zeros(1, 11); % Delta functions timings: delta(t-n*tau) => t = n*tau;
for i=1:11
    delta_func_index(i) = find(t==(i-1)*tau); % t = n*tau;
end
V_tidle = zeros(1, numel(t));
V_tidle(delta_func_index) = 1/delta_t; % V_tidle = 0, except when t = n*tau;
for i=2:numel(t)
   t_i = t(i);
   q_i = V_tidle(i);
   % (V(i+1)-V(i))/step=-V(i)+q(i) => V(i+1)=(-V(i)+q(i))*step+V(i);
   V_euler(i) = V_euler(i-1) - V_euler(i-1)*delta_t/tau + q_i*delta_t/tau;
end
%% Error
Error_s = abs(V_analytical - V_euler);
%% Plot Figure
figure(1)
hold on
yyaxis left
plot(t, V_analytical, 'r', 'linewidth', 2, 'DisplayName', 'V_Analytical')
plot(t, V_euler, 'g', 'linewidth', 2, 'DisplayName', 'V_Euler')
ylabel('V(t)')
yyaxis right
plot(t, Error_s, 'b', 'linewidth', 0.5)
ylabel('Error')
ylim([0 20E-3])
title('V(t) between 0 and 10\tau')
legend({'Analytical Solution','Euler Solution', 'Error'})
xticks(0: tau: 10*tau);
xticklabels({'0','\tau','2\tau','3\tau','4\tau','5\tau','6\tau','7\tau','8\tau','9\tau','10\tau'});
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'q2_A2.png')