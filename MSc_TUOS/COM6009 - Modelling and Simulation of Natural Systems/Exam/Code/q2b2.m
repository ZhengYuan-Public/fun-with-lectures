%% Susceptible-Infected-Recovered (SIR) model with vaccination
% Parameters
N = 10e6; % Population = 1 million
dt = 1e-3; % Time step
t_all = 180; % Total simulation time (days)
times = 0: dt: t_all; % Discretized time
t_v = 30; % Start time for vaccination
% Boundary conditions
S_0 = N; % Initial condition, susceptible individuals
I_0 = 1000; % Initial condition, infected individuals
R_0 = 0; % Initial condition, recovered individuals
% Equation parameters
k = 1/3; % The average period of infectiousness at 3 days
b = 1/2; % Each infected individual makes possible infection contacts every 2 days
alpha = 3000; % 3000 individuals can be vaccinated per day
% Initialize and set boundary condition
S = zeros(1, numel(times)); S(1) = S_0; s = S/N; 
I = zeros(1, numel(times)); I(1) = I_0; i = I/N;
R = zeros(1, numel(times)); R(1) = R_0; r = R/N;
% Make two copies, copy 1 is with vaccination and copy 2 is without vaccination
[S1, S2] = deal(S); [s1, s2] = deal(s);
[I1, I2] = deal(I); [i1, i2] = deal(i);
[R1, R2] = deal(R); [r1, r2] = deal(r);
%% Euler method
% With vaccination at day 10 and alpha individuals can be vaccinated per day
for index=2:numel(times)
    S1(index) = S1(index-1)+dt*(-b*S1(index-1)/N*I1(index-1));
    i1(index) = i1(index-1)+dt*(b*S1(index-1)/N*i1(index-1)-k*i1(index-1));
    I1(index) = i1(index)*N;
    r1(index) = r1(index-1)+dt*(k*i1(index-1)+alpha/N*heaviside(times(index)-t_v));
    R1(index) = r1(index)*N;
end
% Without vaccination
for index=2:numel(times)
    S2(index) = S2(index-1)+dt*(-b*S2(index-1)/N*I2(index-1));
    i2(index) = i2(index-1)+dt*(b*S2(index-1)/N*i2(index-1)-k*i2(index-1));
    I2(index) = i2(index)*N;
    r2(index) = r2(index-1)+dt*(k*i2(index-1));
    R2(index) = r2(index)*N;
end
%% Plot figures
figure(1); hold on;
plot(times, S1, 'blue', 'linewidth', 2);
plot(times, I1, 'red', 'linewidth', 2);
plot(times, R1, 'green', 'linewidth', 2);
plot([t_v, t_v], [0, N], 'r-.',  'linewidth', 2);
legend({'Susceptible  individuals', 'Infected individuals', 'Recovered individuals', 'Start vaccination'});
title({'With vaccination'}); set(gca,'LooseInset',get(gca,'TightInset')); saveas(gcf,'q2b2_1.png');
 
figure(2); hold on;
plot(times, S2, 'blue', 'linewidth', 2);
plot(times, I2, 'red', 'linewidth', 2);
plot(times, R2, 'green', 'linewidth', 2);
legend({'Susceptible  individuals', 'Infected individuals', 'Recovered individuals'});
title({'Without vaccination'}); set(gca,'LooseInset',get(gca,'TightInset')); saveas(gcf,'q2b2_2.png');
