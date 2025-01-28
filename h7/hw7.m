clear all;   % This clears all workspaces
close all;   % This closes all figures 
clc;         % This clears the command window
format long; 

%% Student name, ID number, and homework_number are required
name = 'Colleen Wang';
id = 'A16880246';
hw_num = 7;

% Problem 1
l1 = 2.7;
l2 = 1; % part a,b
l3 = 2.4;
l4 = 3;
lh = 1.8;
lv = 1.5;

% coefficients in Freudenstein's equation
L1 = l1 / l4;
L2 = l1 / l2;
L3 = (l1^2 + l2^2 - l3^2 + l4^2) / (2 * l2 * l4);

% Determine range of rocker angle
theta4_max = acosd(((l2 - l3)^2 - l1^2 - l4^2) / (2 * l1 * l4));
theta4_min = acosd(((l2 + l3)^2 - l1^2 - l4^2) / (2 * l1 * l4));

% For a given theta, solve Freudenstein's equation for theta4
theta2 = linspace(0, 360, 1000);
theta3 = zeros(size(theta2));
theta4 = zeros(size(theta2));

for n = 1:length(theta2)
    Freudenstein = @(theta4) L3 + L2 * cosd(theta4) ...
                        - L1 * cosd(theta2(n)) - cosd(theta2(n) - theta4);

    theta4(n) = fzero(Freudenstein, theta4_max);

    if l2 * cosd(360 - theta2(n)) > l1 - l4 * cosd(180 - theta4(n))
        theta3(n) = 180 - asind((l4 * sind(theta4(n)) - l2 * sind(theta2(n))) / l3);
    else
        theta3(n) = asind((l4 * sind(theta4(n)) - l2 * sind(theta2(n))) / l3);
    end
end

figure(1);
plot(theta2, theta3, 'b', 'LineWidth', 2);
hold on;
plot(theta2, theta4, 'r', 'LineWidth', 2);
xlabel('\theta_2 (degrees)');
ylabel('\theta_3, \theta_4 (degrees)');
title('Evolution of \theta_3 and \theta_4 as a function of \theta_2');
legend('\theta_3', '\theta_4');
grid on;
p1a = 'See figure 1';


% Calculate angular velocities omega3 and omega4
omega2 = 800 * (2 * pi / 60); % rad/s

% Convert degrees to radians for calculation
theta2_rad = deg2rad(theta2);
theta3_rad = deg2rad(theta3);
theta4_rad = deg2rad(theta4);

% Use numerical differentiation to approximate the derivatives
dtheta2_dt = gradient(theta2_rad);
dtheta3_dt = gradient(theta3_rad);
dtheta4_dt = gradient(theta4_rad);

omega3 = (dtheta3_dt ./ dtheta2_dt) * omega2;
omega4 = (dtheta4_dt ./ dtheta2_dt) * omega2;

figure(2);
plot(theta2, omega3, 'b', 'LineWidth', 2);
hold on;
plot(theta2, omega4, 'r', 'LineWidth', 2);
xlabel('\theta_2 (degrees)');
ylabel('\omega_3, \omega_4 (rad/s)');
title('Evolution of \omega_3 and \omega_4 as a function of \theta_2');
legend('\omega_3', '\omega_4');
grid on;
p1b = 'See figure 2';


desired_swing_angle = 28; % rocker 4

% Define the function to calculate the swing angle of theta4 given L2
swing_angle_function = @(L2) calculate_swing_angle(l1, l3, l4, L2);
L2_guess = 1.0;

% Solve for L2
options = optimset('Display','off'); % Suppress display
L2 = fminsearch(@(L2) abs(swing_angle_function(L2) - desired_swing_angle), L2_guess, options);
p1c = L2;


% Calculate the coordinates of point P
xP = l2 * cosd(theta2) + lv * cosd(theta3);
yP = l2 * sind(theta2) + lv * sind(theta3);

% Create the figure
figure(3);
plot(xP, yP, 'b', 'LineWidth', 2);
xlabel('x (m)');
ylabel('y (m)');
title('Trajectory of point P on the coupler');
grid on;
p1d = 'See figure 3';


% Calculate the velocity components (vx and vy) of point P
vx = gradient(xP) ./ gradient(theta2_rad) * omega2;
vy = gradient(yP) ./ gradient(theta2_rad) * omega2;

figure(4);
plot(theta2, vx, 'b', 'LineWidth', 2);
hold on;
plot(theta2, vy, 'r', 'LineWidth', 2);
xlabel('\theta_2 (degrees)');
ylabel('Velocity (m/s)');
title('Evolution of velocity components v_x and v_y as a function of \theta_2');
legend('v_x', 'v_y');
grid on;
p1e = 'See figure 4';


% Function to calculate the swing angle of theta4 given L2
function swing_angle = calculate_swing_angle(l1, l3, l4, L2)
    theta2 = linspace(0, 360, 1000); % 0 to 360 degrees
    L1 = l1 / l4;
    L2_coefficient = l1 / L2;
    L3_coefficient = (l1^2 + L2^2 - l3^2 + l4^2) / (2 * L2 * l4);

    theta4 = zeros(size(theta2));
    for n = 1:length(theta2)
        Freudenstein = @(theta4) L3_coefficient + L2_coefficient * cosd(theta4) ...
                            - L1 * cosd(theta2(n)) - cosd(theta2(n) - theta4);

        theta4(n) = fzero(Freudenstein, 0); % Solve for theta4
    end

    % Calculate the swing angle of theta4
    swing_angle = max(theta4) - min(theta4);
end
