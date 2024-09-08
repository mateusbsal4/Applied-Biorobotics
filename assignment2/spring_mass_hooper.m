%% initialization
clear all;
close all;
clc;

%% Parameters
l_0 = 1;     % [m]
h_max = 1.2; % [m]
g = 9.81;   % [m/s^2]
k = 100;     % [N/m]
m = 1;       %[kg]

flight_mode = 1;        %starts at flight mode, then the check x<1m is performed at the end of every iteration of the for loop

C1 = 0;  %initial condition dx(0) = 0
C2 = h_max; %initial condition x(0) = h_max




timesteps = 0:0.05:10;
x = zeros(size(timesteps));
dx = zeros(size(timesteps));
K = zeros(size(timesteps));
V_grav = zeros(size(timesteps));
V_el = zeros(size(timesteps));
E_total = zeros(size(timesteps));


for i = 1:length(timesteps)
    t = timesteps(i);
    if flight_mode
        dx(i) = -g*t + C1;
        x(i) = -g*(t^2)/2 +C1*t +C2;
        V_el(i) = 0;
    else            %stance
        x(i) = C1*cos(sqrt(k/m)*t)+C2*sin(sqrt(k/m)*t)+l_0-m*g/k;
        dx(i) = -C1*sqrt(k/m)*sin(sqrt(k/m)*t) + C2*sqrt(k/m)*cos(sqrt(k/m)*t);
        V_el(i) = 0.5*k*(x(i)-l_0)^2;
    end
    K(i) = m*dx(i)^2/2;
    V_grav(i) = m*g*x(i);
    E_total(i) = V_el(i) + V_grav(i) + K(i);
    if i>1 && x(i) <= l_0 && x(i-1) > l_0   %flight to stance
        dx_0 = -sqrt(2*g*(h_max-l_0));
        M = [cos(sqrt(k/m)*t) sin(sqrt(k/m)*t); -sqrt(k/m)*sin(sqrt(k/m)*t) sqrt(k/m)*cos(sqrt(k/m)*t)];
        v = [m*g/k; dx_0];
        C = M \v;
        C1 = C(1);
        C2 = C(2);
        flight_mode = 0;

    elseif i> 1 && x(i) >l_0 && x(i-1) <=l_0      %stance to flight 
        x_0 = l_0;
        dx_0 = sqrt(2*g*(h_max-l_0));
        C1 = dx_0 +g*t;
        C2 = x_0 +g*t^2/2 -C1*t;
        flight_mode = 1;


    end
end

clf
figure(1)
subplot(2, 1, 1)
plot(timesteps, x, 'LineWidth', 2, 'Color', 'b')
grid on
xlabel('t [s]', 'FontSize', 14)
ylabel('x [m]', 'FontSize', 14)
title("Displacement of the mass over time", 'FontSize', 15)
subplot(2, 1, 2)
plot(timesteps, dx, 'LineWidth', 2, 'Color', 'r')
grid on
xlabel('t [s]', 'FontSize', 14')
ylabel('$\dot{x} [m/s]$','Interpreter','Latex', 'FontSize', 14)
title("Velocity of the mass over time", 'FontSize', 15)


figure(2)
plot(timesteps, K, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Kinetic Energy')
hold on
plot(timesteps, V_el, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Elastic Potential Energy')
plot(timesteps, V_grav, 'LineWidth', 2, 'Color', 'g', 'DisplayName', 'Gravitational Potential Energy')
plot(timesteps, E_total, 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Total Energy')
grid on
hold off
xlabel('t [s]', 'FontSize', 14)
ylabel('Energy [J]', 'FontSize', 14)
title("Contributions of the gravitational, elastic and kinetic energies of the system over time", 'FontSize', 15)
legend('FontSize', 11)


