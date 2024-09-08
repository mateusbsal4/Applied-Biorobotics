clear;
close all;
clc;

function K = F(X, stance_1, stance_2, l0, g, m, k, x_FP1, x_FP2)
    K = zeros(size(X));     % both 1x4 vectors
    x = X(1);
    dx = X(2);
    y = X(3);
    dy = X(4);
    assert(stance_1 + stance_2 > 0, "Simulation failure: Model does not account for running!")   
    K(1) = dx;
    if stance_1 && stance_2  % Both legs in stance
        K(2) = -(k/m) * ((x_FP1-x) * (l0 / sqrt((x-x_FP1)^2 + y^2) - 1) + (x_FP2-x) * (l0 / sqrt((x-x_FP2)^2 + y^2) - 1));
    elseif stance_1 && ~stance_2  % Only 1st leg in stance
        K(2) = -(k/m) * (x_FP1-x) * (l0 / sqrt((x-x_FP1)^2 + y^2) - 1);
    elseif ~stance_1 && stance_2  % Only 2nd leg in stance
        K(2) = -(k/m) * (x_FP2-x) * (l0 / sqrt((x-x_FP2)^2 + y^2) - 1);
    end
    K(3) = dy;
    if stance_1 && stance_2  
        K(4) = -g + (k/m) * y * (l0 / sqrt((x-x_FP1)^2 + y^2) + l0 / sqrt((x-x_FP2)^2 + y^2) - 2);
    elseif stance_1 && ~stance_2
        K(4) = -g + (k/m) * y * (l0 / sqrt((x-x_FP1)^2 + y^2) - 1);
    elseif ~stance_1 && stance_2
        K(4) = -g + (k/m) * y * (l0 / sqrt((x-x_FP2)^2 + y^2) - 1);
    end
end

function [x1, x2, x3, x4] = rk4_solver(X, stance_1, stance_2, l0, g, m, k, x_FP1, x_FP2, h)  %%implements rk4
    K1 = F(X, stance_1, stance_2, l0, g, m, k, x_FP1, x_FP2);
    K2 = F(X + 0.5 * h * K1, stance_1, stance_2, l0, g, m, k, x_FP1, x_FP2);
    K3 = F(X + 0.5 * h * K2, stance_1, stance_2, l0, g, m, k, x_FP1, x_FP2);
    K4 = F(X + h * K3, stance_1, stance_2, l0, g, m, k, x_FP1, x_FP2);
    X = X + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4);
    x1 = X(1);
    x2 = X(2);
    x3 = X(3);
    x4 = X(4);    
end

%% Constants definition
g = 9.81;
k = 10000;
m = 80;
alpha_TD = 68 * pi / 180;
l0 = 1;

% Initialization x0, dx0, y0, dy0
x = 0;
y = 0.97;
dx = 1.05;
dy = 0; 
l1 = y;
l2 = l0;
x_FP1 = 0;
x_FP2 = 0.243; % to ensure that l2 > l0 initially
stance_1 = 1;
stance_2 = 0;
first_touchdown = 0;
h = 0.001; % Time step
timesteps = 0:h:10;
% Vectors containing the x y displacements, forces and energies at each timestep
x_vec = zeros(size(timesteps));
y_vec = zeros(size(timesteps)); 
dy_vec = zeros(size(timesteps));
dx_vec = zeros(size(timesteps));
f_vert1 = zeros(size(timesteps));
f_horiz1 = zeros(size(timesteps));
f_vert2 = zeros(size(timesteps));
f_horiz2 = zeros(size(timesteps));
K = zeros(size(timesteps));
V_el = zeros(size(timesteps));
V_grav = zeros(size(timesteps));
E_total = zeros(size(timesteps));
% Touchdown instants for spring 1
td1_ts = [];
%Touchdown instants for spring 2
td2_ts = [];
%Takeoff  instants for spring 1  
to1_ts = [];
%Takeoff instants for spring 2  
to2_ts = [];

for i = 1:length(timesteps)
    if y <= l0 * sin(alpha_TD) && ~stance_1 % Touchdown of spring 1
        x_FP1 = x + l0 * cos(alpha_TD);
        stance_1 = 1;
    elseif y <= l0 * sin(alpha_TD) && ~stance_2  % Touchdown of spring 2
        x_FP2 = x + l0 * cos(alpha_TD);
        first_touchdown = 1;
        stance_2 = 1;
    end
    if stance_1
        l1 = sqrt((x - x_FP1) ^ 2 + y ^ 2);
    else
        l1 = l0;
    end
    if stance_2
        l2 =  sqrt((x - x_FP2) ^ 2 + y ^ 2);
    else 
        l2 = l0;
    end
    if l1 >=l0 
        if stance_1
            if isempty(to1_ts)         
                to1_ts = [to1_ts, timesteps(i)];
            elseif timesteps(i)-td1_ts(length(td1_ts)) > 0.09 && timesteps(i)-to1_ts(length(to1_ts)) > 0.09       
                to1_ts = [to1_ts, timesteps(i)];
            end
        end
        stance_1 = 0;
        l1 = l0;
        f_horiz1_i = 0;
        f_vert1_i = 0;
        %E_el_1 = 0;
    else
        if x_FP1 ==x + l0 * cos(alpha_TD) 
            if isempty(td1_ts)         
                td1_ts = [td1_ts, timesteps(i)]; 
            elseif  timesteps(i)-to1_ts(length(to1_ts)) > 0.09  && timesteps(i)-td1_ts(length(td1_ts)) > 0.09         
                td1_ts = [td1_ts, timesteps(i)];                    
            end   
        end
        stance_1 = 1;
        f_horiz1_i = k*(l0-l1)*((x-x_FP1)/l1);
        f_vert1_i = k*(l0-l1)*(y/l1);
        %E_el_1 = 0.5*k*(l1-l0)^2;
    end 
    if l2 >= l0
        if stance_2 
            if isempty(to2_ts)         
                to2_ts = [to2_ts, timesteps(i)];
            elseif timesteps(i)-td2_ts(length(td2_ts)) > 0.09 && timesteps(i)-to2_ts(length(to2_ts)) > 0.09        
                to2_ts = [to2_ts, timesteps(i)];
            end
        end
        stance_2 = 0;
        l2 = l0;
        f_horiz2_i = 0;
        f_vert2_i = 0;
        %E_el_2 = 0;
    elseif l2<=l0 && first_touchdown
        if x_FP2 ==x + l0 * cos(alpha_TD) 
            if isempty(td2_ts)          
                td2_ts = [td2_ts, timesteps(i)];
            elseif timesteps(i)-to2_ts(length(to2_ts)) > 0.09 && timesteps(i)-td2_ts(length(td2_ts)) > 0.09      
                td2_ts = [td2_ts, timesteps(i)];
            end
        end
        stance_2 = 1;
        f_horiz2_i = k*(l0-l2)*((x-x_FP2)/l2);
        f_vert2_i = k*(l0-l2)*(y/l2);
        %E_el_2 = 0.5*k*(l2-l0)^2;
    elseif l2<=l0 && ~first_touchdown
        stance_2 = 0;
        f_horiz2_i = k*(l0-l2)*((x-x_FP2)/l2);
        f_vert2_i = k*(l0-l2)*(y/l2);
        %E_el_2 = 0.5*k*(l2-l0)^2;
    end
    X = [x; dx; y; dy]; % Ensure X is a column vector
    [x, dx, y, dy] = rk4_solver(X, stance_1, stance_2, l0, g, m, k, x_FP1, x_FP2, h);
    x_vec(i) = x;
    y_vec(i) = y;
    dx_vec(i) = dx;
    dy_vec(i) = dy;
    f_vert1(i) = f_vert1_i;
    f_horiz1(i) = f_horiz1_i;
    f_vert2(i) = f_vert2_i;
    f_horiz2(i) = f_horiz2_i;    
    K(i) = 0.5*m*(dx^2+dy^2);
    E_el_2 = 0.5*k*(l2-l0)^2;
    E_el_1 = 0.5*k*(l1-l0)^2;
    V_el(i) = E_el_1 + E_el_2;
    V_grav(i) = m*g*y;
    %E_total_i = E_el_1 + E_el_2 + m*g*y + 0.5*m*(dx^2+dy^2);
    E_total(i) = K(i) +V_el(i) +V_grav(i);
end


figure(1)
plot(x_vec, y_vec, 'LineWidth', 2, 'Color', 'b')
xlabel('x [m]', 'FontSize', 16)
ylabel('y [m]', 'FontSize', 16)
title("Trajectory of the mass",  'FontSize', 17)




figure(2)
subplot(2, 1, 1)
plot(timesteps, f_horiz1, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Left leg')
hold on
plot(timesteps, f_horiz2, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Right leg')
grid on
hold on
scatter(td1_ts, zeros(size(td1_ts)), 'filled', 'Marker', '^', 'MarkerFaceColor', 'b', 'DisplayName', 'Left TO')                %left leg takeoff and touchdown markers
hold on
scatter(to1_ts, zeros(size(to1_ts)), 'filled', 'Marker', 'v', 'MarkerFaceColor', 'b', 'DisplayName', 'Left TD')
hold on
scatter(td2_ts, zeros(size(td2_ts)), 'filled', 'Marker', '^', 'MarkerFaceColor', 'k', 'DisplayName', 'Right TO')                 %right leg takeoff and touchdown markers
hold on
scatter(to2_ts, zeros(size(to2_ts)), 'filled', 'Marker', 'v', 'MarkerFaceColor', 'k', 'DisplayName', 'Right TD')
hold off
xlabel('t [s]', 'FontSize', 16)
ylabel('$F_{el, x} [N]$','Interpreter','Latex', 'FontSize', 16)
title("Horizontal elastic forces over time", 'FontSize', 17)
legend('FontSize', 11)
subplot(2, 1, 2)
plot(timesteps, f_vert1, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Left leg')
hold on
plot(timesteps, f_vert2, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Right leg')
grid on
hold off
xlabel('t [s]', 'FontSize', 16)
ylabel('$F_{el, y} [N]$', 'Interpreter', 'Latex', 'FontSize', 16)
title("Vertical elastic forces over time", 'FontSize', 17)
legend('FontSize', 11)

figure(3)
plot(y_vec,dy_vec , 'LineWidth', 2, 'Color', 'b')
xlabel('y [m]', 'FontSize', 16)
ylabel('$\dot{y}$ [m/s]', 'Interpreter', 'Latex', 'FontSize', 16)  % Corrected label unit to N
title("Phase plot", 'FontSize', 17)


figure(4)
plot(timesteps, K, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Kinetic Energy')
hold on
plot(timesteps, V_el, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Elastic Potential Energy')
plot(timesteps, V_grav, 'LineWidth', 2, 'Color', 'g', 'DisplayName', 'Gravitational Potential Energy')
plot(timesteps, E_total, 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Total Energy')
grid on
hold off
xlabel('t [s]', 'FontSize', 16)
ylabel('Energy [J]', 'FontSize', 16)
title("Contributions of the gravitational, elastic and kinetic energies of the system over time", 'FontSize', 17)
legend('FontSize', 11)


%Average velocity computation
window_size = length(timesteps)/20; 
moving_avg_dx = movmean(dx_vec, window_size);
fprintf('Maximum instantaneous forward speed: %.3f\n', max(dx_vec));
fprintf('Maximum average forward speed: %.3f\n', max(moving_avg_dx));
figure(5)
plot(timesteps,dx_vec , 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Instantaneous velocity')
hold on;
plot(timesteps, moving_avg_dx, 'r', 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Average velocity');  % Moving average
hold off;
xlabel('t [s]', 'FontSize', 16)
ylabel('$\dot{x}$ [m/s]', 'Interpreter', 'Latex', 'FontSize', 16)  
title("Forward velocity of the mass over time", 'FontSize', 17)
legend('FontSize', 11)