clear
close all
%% Boundary conditions
t_start = 0; t_end = 1000;
v_start = -3; w_start = -3;
yvec0 = [v_start w_start];
%% Solve with ode45
options = odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5]);
[tvec, ymat] = ode45(@fitzhugh_nagumo_model, [t_start t_end], yvec0, options);
%% Plot figure
figure(1)
hold on;
plot(ymat(:,1), ymat(:,2),'LineWidth', 3)
ylabel('W(t)');
xlabel('V(t)');
title('Phase plane plot')
axis tight

N=30; % number of samples
xmin=-4; xmax=4;  
ymin=-4; ymax=4;

[xmat,ymat] = meshgrid(xmin: (xmax-xmin)/N: xmax, ymin: (ymax-ymin)/N: ymax); 
dxmat = zeros(size(xmat));
dymat = zeros(size(ymat));
for xi=1:N
    for yi=1:N
        yprime = fitzhugh_nagumo_model(0, [xmat(xi,yi);ymat(xi,yi)]);
        dxmat(xi,yi) = yprime(1,1);
        dymat(xi,yi) = yprime(2,1);
    end
end
arrowsplot=quiver(xmat, ymat, dxmat, dymat); 
set(arrowsplot, 'linewidth', 1.5);
arrowsplot.AutoScaleFactor = 1.5;
% Plot nullclines
I = 0.5; e = 1/12.5; b0 = 0.8; b1 = 0.7;

syms x  y
V_nullclines = x - x.^3/3 + I;
W_nullclines = b1*x + b0;
p1 = fplot(V_nullclines, 'r', 'linewidth', 2);
p2 = fplot(W_nullclines, 'g', 'linewidth', 2);

% Find the intercestion of V_nullclines and W_nullclines
[X, Y] = solve(V_nullclines == y, W_nullclines == y);
X = double(X); % Intersection x coordinates
Y = double(Y); % Intersection y coordinates
p3 = plot(X(1), Y(1), 'k.', 'MarkerSize', 20);
legend([p1(1), p2(1), p3(1)], 'V-nullcline', 'W-nullcline', 'Equilibrium point')
xlim([-4 4]);
ylim([-4 4]);
set(gca,'LooseInset',get(gca,'TightInset'));
saveas(gcf,'q2_B_1.png')
%% Fitzhugh-Nagumo Model
function yprime = fitzhugh_nagumo_model(~, yvec)
    % yvec = [v0, w0]
    % yprime = [dvdt, dwdt]
    % Parameters
    I = 0.5; e = 1/12.5; b0 = 0.8; b1 = 0.7;
    yprime(1, 1) = yvec(1) - yvec(1).^3/3 - yvec(2) + I;
    yprime(2, 1) = e.*(b0 + b1.*yvec(1) - yvec(2));
end