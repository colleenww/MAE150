clear all;   % This clears all workspaces
close all;   % This closes all figures 
clc;         % This clears the command window
format long; 

%% Student name, ID number, and homework_number are required
name = 'Colleen Wang';
id = 'A16880246';
hw_num = 3;

%% Problem 1
x = linspace(-2, 2, 100);
y = x.^2;
offsets = [0.25 0.5 0.75 1];

figure(1); hold on; grid on;
plot(x, y, 'k-', 'LineWidth', 2);
for i = 1:length(offsets)
    D = offsets(i);
    x_offsets = x - 2*x./(sqrt(4*x.^2 + 1)) * D;
    y_offsets = y + 1./(sqrt(4*x.^2 + 1)) * D;
    plot(x_offsets, y_offsets, 'LineWidth', 1.5);
end
xlabel('x'); ylabel('y');
title('Parabola with Offset Curves');
legend('y = x^2','D = 0.25','D = 0.5','D = 0.75','D=1','Location','northeast');
p1a = 'See figure 1';


dy = 2*x; dy2 = 2;
K = dy2 ./ ((1 + dy.^2).^(3/2));
R = 1 ./ K;

figure(2); hold on;
subplot(2, 1, 1);
plot(x, R, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('Radius of Curvature (R)');
legend('R');
title('Radius of Curvature of y = x^2');

subplot(2, 1, 2);
plot(x, K, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('Curvature (K)');
legend('K');
title('Curvature of y = x^2');

p1b = 'See figure 2';

p1c = 'When the offset distance D exceeds the radius of curvature (R), the offset curves can self-intersect and create loops with sharp points. This can be seen by the curves of y=x^2 offset by a distance of D = 0.75 and D = 1 because as the curve moves "upwards", the radius of curvature is decreasing.';



%% Problem 2
x = {-1, 0, 0, 1, 1, -1};
y = {0, 1, 1, 0, 0, 0};
dydx = {10, -10, -10, -10, 5, -5};

% cubic hermite spline interpolation
figure (3); hold on;
for i = 1:3
    xi = [x{2*i-1}, x{2*i}];
    yi = [y{2*i-1}, y{2*i}];
    dydxi = [dydx{2*i-1}, dydx{2*i}];
    
    t = linspace(0, 1, 100);
    h00 = 2*t.^3 - 3*t.^2 + 1;
    h10 = t.^3 - 2*t.^2 + t;
    h01 = -2*t.^3 + 3*t.^2;
    h11 = t.^3 - t.^2;
    f = h00*yi(1) + h10*(xi(2)-xi(1))*dydxi(1) + h01*yi(2) + h11*(xi(2)-xi(1))*dydxi(2);
    
    plot(linspace(xi(1), xi(2), 100), f, 'LineWidth', 2);
end

% Formatting the Plot
title('Cubic Hermite Spline Interpolation');
xlabel('x');
ylabel('y');
legend('curve #1', 'curve #2', 'curve #3');
grid on;

p2 = 'See figure 3';



%% Problem 3
% Cam parameters
r_follower = 6e-3;      % [m] follower radius
r_prime = 0.02;         % [m] radius of prime circle
omega = 1000*2*pi/60;   % [rad/sec] convert CAM rotational speed from RPM to rad/sec
r_base = r_prime - r_follower;
dt = deg2rad(2);        % radial resolution [rad];

% Define various displacement profiles using symbolic math
syms theta beta L
dwell = 0;
polynomial = L*(10/beta^3*theta^3 - 15/beta^4*theta^4 + 6/beta^5*theta^5);
cycloidal = L*(theta/beta - 1/(2*pi)*sin(2*pi*theta/beta));
harmonic = (L/2)*(1 - cos((pi*theta)/beta));

% Note:
% v = dy/dt = dy/dtheta * omega
% a = d2y/dt2= d2y/dtheta2 * omega^2
% j = d3y/dt3 = d2y/dtheta3 * omega^3

% Generate cam surface for each section
% 3-4-5 polynomial rise (0 < theta < 110) 0mm to 10mm
L1 = 0.01; 
beta1 = deg2rad(110);
t1 = 0:dt:beta1;
y1 = subs(polynomial,{theta,beta,L},{t1,beta1,L1});
v1 = subs(diff(polynomial,theta)*omega,{theta,beta,L},{t1,beta1,L1});
a1 = subs(diff(polynomial,theta,2)*omega^2,{theta,beta,L},{t1,beta1,L1});
j1 = subs(diff(polynomial,theta,3)*omega^3,{theta,beta,L},{t1,beta1,L1});

% Dwell (110 < theta < 120)
L2 = L1; % constant displacement
beta2 = deg2rad(120);
t2 = beta1:dt:beta2;
y2 = y1(end) + subs(dwell,{theta,beta,L},{t2-t2(1),beta2-beta1,L2-L1});
v2 = subs(diff(dwell,theta)*omega,{theta,beta,L},{t2-t2(1),beta2-beta1,L2-L1});
a2 = subs(diff(dwell,theta,2)*omega^2,{theta,beta,L},{t2-t2(1),beta2-beta1,L2-L1});
j2 = subs(diff(dwell,theta,3)*omega^3,{theta,beta,L},{t2-t2(1),beta2-beta1,L2-L1});

% Cycloidal rise (120 < theta < 200) 10 mm to 20mm
L3 = 0.02;
beta3 = deg2rad(200);
t3 = beta2:dt:beta3;
y3 = y2(end) + subs(cycloidal,{theta,beta,L},{t3-t3(1),beta3-beta2,L3-L2});
v3 = subs(diff(cycloidal,theta)*omega,{theta,beta,L},{t3-t3(1),beta3-beta2,L3-L2});
a3 = subs(diff(cycloidal,theta,2)*omega^2,{theta,beta,L},{t3-t3(1),beta3-beta2,L3-L2});
j3 = subs(diff(cycloidal,theta,3)*omega^3,{theta,beta,L},{t3-t3(1),beta3-beta2,L3-L2});

% Dwell (200 < theta < 220)
L4 = L3; % constant displacement
beta4 = deg2rad(220);
t4 = beta3:dt:beta4;
y4 = y3(end) + subs(dwell,{theta,beta,L},{t4-t4(1),beta4-beta3,L4-L3});
v4 = subs(diff(dwell,theta)*omega,{theta,beta,L},{t4-t4(1),beta4-beta3,L4-L3});
a4 = subs(diff(dwell,theta,2)*omega^2,{theta,beta,L},{t4-t4(1),beta4-beta3,L4-L3});
j4 = subs(diff(dwell,theta,3)*omega^3,{theta,beta,L},{t4-t4(1),beta4-beta3,L4-L3});

% Harmonic (220 < theta < 340) 20mm to 0mm
L5 = 0;
beta5 = deg2rad(340);
t5 = beta4:dt:beta5;
y5 = y4(end) + subs(harmonic,{theta,beta,L},{t5-t5(1),beta5-beta4,L5-L4});
v5 = subs(diff(harmonic,theta)*omega,{theta,beta,L},{t5-t5(1),beta5-beta4,L5-L4});
a5 = subs(diff(harmonic,theta,2)*omega^2,{theta,beta,L},{t5-t5(1),beta5-beta4,L5-L4});
j5 = subs(diff(harmonic,theta,3)*omega^3,{theta,beta,L},{t5-t5(1),beta5-beta4,L5-L4});

% Dwell (340 < theta < 360)
L6 = L5; % constant displacement
beta6 = deg2rad(360);
t6 = beta5:dt:beta6;
y6 = y5(end) + subs(dwell,{theta,beta,L},{t6-t6(1),beta6-beta5,L6-L5});
v6 = subs(diff(dwell,theta)*omega,{theta,beta,L},{t6-t6(1),beta6-beta5,L6-L5});
a6 = subs(diff(dwell,theta,2)*omega^2,{theta,beta,L},{t6-t6(1),beta6-beta5,L6-L5});
j6 = subs(diff(dwell,theta,3)*omega^3,{theta,beta,L},{t6-t6(1),beta6-beta5,L6-L5});

% Combine all sections
t = [t1, t2, t3, t4, t5, t6];
dis = eval([y1, y2, y3, y4, y5, y6]);       % eval: convert symbolic type to double precision
vel = eval([v1, v2, v3, v4, v5, v6]);
acc = eval([a1, a2, a3, a4, a5, a6]);
jerk = eval([j1, j2, j3, j4, j5, j6]);

figure(4); hold on;    
subplot(4,1,1);
plot(rad2deg(t), dis*1000,'k');         % displacement
title('Displacement Profile');
xlabel('\theta [degrees]');
ylabel('displacement [mm]');
box on; grid on;
set(gca,'XTick',0:30:360); hold off;

subplot(4,1,2);
plot(rad2deg(t), vel*1000, 'k');        % velocity
title('Velocity Profile');
xlabel('\theta [degrees]');
ylabel('velocity [mm/s]');
box on; grid on;
set(gca,'XTick',0:30:360); hold off;

subplot(4,1,3);
plot(rad2deg(t), acc*1000, 'k');        % acceleration
title('Acceleration Profile');
xlabel('\theta [degrees]');
ylabel('acceleration [mm/s^2]'); 
box on; grid on;
set(gca,'XTick',0:30:360); hold off;

subplot(4,1,4);
plot(rad2deg(t), jerk*1000, 'k');       % jerk
title('Jerk Profile');
xlabel('\theta [degrees]');
ylabel('jerk [mm/s^3]'); 
box on; grid on;
set(gca,'XTick',0:30:360); hold off;

p3a = 'See figure 4';


% Define parameters
r_follower = 6e-3;        % [m] Follower radius
r_prime = 0.02;           % [m] Radius of prime circle
omega = 1000*2*pi/60;     % [rad/sec] Convert CAM rotational speed from RPM to rad/sec
r_base = r_prime - r_follower;  % [m] Base circle radius
dt = deg2rad(2);          % [rad] Radial resolution

% Angular domain
t = 0:dt:2*pi;            % Complete rotation in radians

% Initialize displacement array
dis = zeros(size(t));

% 3-4-5 Polynomial rise (0 to 110 degrees)
indices = t <= deg2rad(110);
dis(indices) = 10e-3 * (3 * (t(indices)/deg2rad(110)).^2 - 2 * (t(indices)/deg2rad(110)).^3);

% Dwell (110 to 120 degrees)
indices = t > deg2rad(110) & t <= deg2rad(120);
dis(indices) = dis(t == deg2rad(110));

% Cycloidal rise (120 to 200 degrees)
indices = t > deg2rad(120) & t <= deg2rad(200);
t_normalized = (t(indices) - deg2rad(120)) / (deg2rad(200) - deg2rad(120));
dis(indices) = 10e-3 + 10e-3 * (t_normalized - (sin(2*pi*t_normalized)/(2*pi)));

% Dwell (200 to 220 degrees)
indices = t > deg2rad(200) & t <= deg2rad(220);
dis(indices) = dis(t == deg2rad(200));

% Harmonic fall (220 to 340 degrees)
indices = t > deg2rad(220) & t <= deg2rad(340);
t_normalized = (t(indices) - deg2rad(220)) / (deg2rad(340) - deg2rad(220));
dis(indices) = 20e-3 * (1 - 0.5 * (1 - cos(pi * t_normalized)));

% Dwell (340 to 360 degrees)
indices = t > deg2rad(340) & t <= 2*pi;
dis(indices) = 0;

% Transform displacements to Cartesian coordinates
x_pitch = (r_prime + dis) .* sin(t);
y_pitch = (r_prime + dis) .* cos(t);

% Cam contour (accounting for the follower radius)
x_cam = (r_prime + dis + r_follower) .* sin(t);
y_cam = (r_prime + dis + r_follower) .* cos(t);

% Prime circle
x_prime = r_prime * sin(t);
y_prime = r_prime * cos(t);

% Base circle
x_base = r_base * sin(t);
y_base = r_base * cos(t);

% Plotting
figure(5); hold on;
plot(x_prime*1000, y_prime*1000, 'b-', 'DisplayName', 'Prime Circle','LineWidth',2);
plot(x_base*1000, y_base*1000, 'b--', 'DisplayName', 'Base Circle','LineWidth',2);
plot(x_pitch*1000, y_pitch*1000, 'r-', 'DisplayName', 'Pitch Curve','LineWidth',2);
plot(x_cam*1000, y_cam*1000, 'g-', 'DisplayName', 'Cam Contour','LineWidth',2);
legend show; axis equal; grid on;
xlabel('X Coordinate (mm)');
ylabel('Y Coordinate (mm)');
title('Cam Design Profile');

p3b = 'See figure 5';

% Define parameters
r_follower = 5e-3;        % [m] Follower radius
r_prime = 0.02;           % [m] Radius of prime circle
omega = 1000*2*pi/60;     % [rad/sec] Convert CAM rotational speed from RPM to rad/sec
dt = deg2rad(2);          % [rad] Radial resolution
t = 0:dt:2*pi;            % Complete rotation in radians

% Displacement profile initialization
dis = zeros(size(t));

% 3-4-5 Polynomial rise (0 to 110 degrees)
indices = t <= deg2rad(110);
dis(indices) = 10e-3 * (3 * (t(indices)/deg2rad(110)).^2 - 2 * (t(indices)/deg2rad(110)).^3);

% Dwell from 110 to 120 degrees and from 340 to 360 degrees
indices110_120 = t > deg2rad(110) & t <= deg2rad(120);
indices340_360 = t > deg2rad(340) & t <= 2*pi;
dis(indices110_120 | indices340_360) = dis(t == deg2rad(110));

% Cycloidal rise from 120 to 200 degrees
indices = t > deg2rad(120) & t <= deg2rad(200);
theta_rel = t(indices) - deg2rad(120);
T = deg2rad(200) - deg2rad(120);
dis(indices) = 10e-3 + 10e-3 * (theta_rel/T - sin(2*pi*theta_rel/T)/(2*pi));

% Dwell from 200 to 220 degrees
indices200_220 = t > deg2rad(200) & t <= deg2rad(220);
dis(indices200_220) = dis(t == deg2rad(200));

% Harmonic fall from 220 to 340 degrees
indices = t > deg2rad(220) & t <= deg2rad(340);
t_normalized = (t(indices) - deg2rad(220)) / (deg2rad(340) - deg2rad(220));
dis(indices) = 20e-3 * (1 - 0.5 * (1 - cos(pi * t_normalized)));

% Compute velocity using finite differences
vel = [0, diff(dis)] / dt;  % Simple approximation of derivative

% Calculate pressure angle phi
denom = sqrt(r_prime^2) .* dis;
denom(denom == 0) = 1e-6;  % Prevent division by zero
phi = atan((vel / omega) ./ denom);

% Plotting
figure(6); clf; hold on;
plot(rad2deg(t), rad2deg(phi), 'k');       % Pressure angle
tmp = [rad2deg(t) fliplr(rad2deg(t))];
patch(tmp, [-30*ones(1,numel(t)) 30*ones(1,numel(t))], 'b','FaceAlpha', 0.2);
ylabel('\phi [degrees]'); 
xlabel('\theta [degrees]');
axis tight; box on; grid on;
set(gca, 'XTick', 0:30:360, 'XTickLabel', 0:30:360);
hold off;

p3c = 'See figure 6';
p3d = 'This is not a good CAM design because of the acceleration and jerk discontinuities. It can be amended by modifying the cam profile to smooth transitions between phases.';
