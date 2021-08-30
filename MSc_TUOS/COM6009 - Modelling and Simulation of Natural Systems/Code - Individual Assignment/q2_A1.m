clear
close all
%% Parameters
tau = 1;
period = 2*tau;
num_periods = 5; % Number of periods;
delta_t = 0.01; % Step size
t = (0: delta_t: num_periods*period); % Discretized time;
omega = pi/tau; % The period of square() is T=2*pi/omega=2*tau => omega = pi/tau;
%% Analytical Method
% Boundary condition: tau = 1; V(0) = 0;
V_analytical= zeros(1, numel(t));
f = @(s, t) exp(s-t).*square(omega*s);
for i=1:numel(t)
    t_i = t(i);
    V_analytical(i) = integral(@(s) f(s, t_i), 0, t_i);
end
%% Euler's method
V_tidle = square(omega*t);
V_euler = zeros(1, numel(t));
V_euler(1) = 0; % Same boundary condition as the analytical method;
for i=2:numel(t)
   t_i = t(i);
   % (V(i+1)-V(i))/step=-V(i)+q(i) => V(i+1)=(-V(i)+q(i))*step+V(i);
   V_euler(i) = (-V_euler(i-1)+V_tidle(i-1))*delta_t + V_euler(i-1);
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
saveas(gcf,'q2_A1.png')