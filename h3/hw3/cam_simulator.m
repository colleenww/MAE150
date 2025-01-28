clear all; close all; clc;

%% Simulation parameters:
export_cam_contour = false;
cam_contour_filename = 'cam1_profile.txt';
export_animation = false;
animation_filename = 'cam1_animation.gif';

%% Cam parameters
r_follower = 6e-3;      % [m] follower radius
r_prime = 0.02;         % [m] radius of prime circle
omega = 1000*2*pi/60;   % [rad/sec] convert CAM rotational speed from RPM to rad/sec
r_base = r_prime - r_follower;
dt = deg2rad(2);        % radial resolution [rad];

%% Define various displacement profiles using symbolic math
syms theta beta L
dwell = 0;
polynomial = L*(10/beta^3*theta^3 - 15/beta^4*theta^4 + 6/beta^5*theta^5);
cycloidal = L*(theta/beta - 1/(2*pi)*sin(2*pi*theta/beta));

% Note:
% v = dy/dt = dy/dtheta * omega
% a = d2y/dt2= d2y/dtheta2 * omega^2
% j = d3y/dt3 = d2y/dtheta3 * omega^3

%% Generate cam surface for each section
% Cycloidal rise (0 < theta < 120) 0 mm to 10mm
L1 = 0.01;
beta1 = deg2rad(120);
t1 = 0:dt:beta1-dt;
y1 = subs(cycloidal,{theta,beta,L},{t1,beta1,L1});
v1 = subs(diff(cycloidal,theta)*omega,{theta,beta,L},{t1,beta1,L1});
a1 = subs(diff(cycloidal,theta,2)*omega^2,{theta,beta,L},{t1,beta1,L1});
j1 = subs(diff(cycloidal,theta,3)*omega^3,{theta,beta,L},{t1,beta1,L1});

% Dwell (120 < 0 < 300)
L2 = 0.01;
beta2 = deg2rad(300);
t2 = beta1:dt:beta2-dt;
y2 = y1(end) + subs(dwell,{theta,beta,L},{t2-t2(1),beta2-beta1,L2-L1});
v2 = subs(diff(dwell,theta)*omega,{theta,beta,L},{t2-t2(1),beta2-beta1,L2-L1});
a2 = subs(diff(dwell,theta,2)*omega^2,{theta,beta,L},{t2-t2(1),beta2-beta1,L2-L1});
j2 = subs(diff(dwell,theta,3)*omega^3,{theta,beta,L},{t2-t2(1),beta2-beta1,L2-L1});

% 3-4-5 polynomial fall (300 < 0 < 360) 10mm to 0mm
L3 = 0;
beta3 = deg2rad(360);
t3 = beta2:dt:beta3-dt;
y3 = y2(end) + subs(polynomial,{theta,beta,L},{t3-t3(1),beta3-beta2,L3-L2});
v3 = subs(diff(polynomial,theta)*omega,{theta,beta,L},{t3-t3(1),beta3-beta2,L3-L2});
a3 = subs(diff(polynomial,theta,2)*omega^2,{theta,beta,L},{t3-t3(1),beta3-beta2,L3-L2});
j3 = subs(diff(polynomial,theta,3)*omega^3,{theta,beta,L},{t3-t3(1),beta3-beta2,L3-L2});

% Combine all sections
t = [t1, t2, t3];
dis = eval([y1, y2, y3]);       % eval: convert symbolic type to double precision
vel = eval([v1, v2, v3]);
acc = eval([a1, a2, a3]);
jerk = eval([j1, j2, j3]);

%% Construct cam geometry
% Obtain pressure angle
e = 0;
phi = atan((vel/omega-e)./(sqrt(r_prime^2-e^2)*dis));

% Pitch curve
x_pitch = (r_prime + dis).*sin(t);
y_pitch = (r_prime + dis).*cos(t);

% Cam contour
x_cam = x_pitch - r_follower*sin(t - phi);
y_cam = y_pitch - r_follower*cos(t - phi);

% Base circle
x_base = r_base*sin(t);
y_base = r_base*cos(t);

% Prime circle
x_prime = r_prime*sin(t);
y_prime = r_prime*cos(t);

%% Plot motion profiles and animate cam motion
cam_in(1,:) = x_cam;
cam_in(2,:) = y_cam;
pitch_in(1,:) = x_pitch;
pitch_in(2,:) = y_pitch;

cam_in(:,end) = cam_in(:,1);        % enforce periodicity
pitch_in(:,end) = pitch_in(:,1);

if export_cam_contour
    csvwrite(cam_contour_filename,1000*[cam_in; zeros(1:length(t))]');
end

figure('Units','inches','Position',[1 1 15 10]);
for nt = 1:1:length(t)
    sp(1) = subplot(6,1,1);
    plot(rad2deg(t), dis*1000,'k');         % displacement
    hold on;
    plot(rad2deg(t(nt)), dis(nt)*1000,'ko','MarkerFaceColor','k');
    ylabel('displacement [mm]');
    box on; grid on; axis tight;
    set(gca,'XTick',0:30:360,'XTickLabel',[]); hold off;

    sp(2) = subplot(6,1,2);
    plot(rad2deg(t), vel*1000, 'k');        % velocity
    hold on;
    plot(rad2deg(t(nt)), vel(nt)*1000,'ko','MarkerFaceColor','k');
    ylabel('velocity [mm/s]'); axis tight;
    box on; grid on;
    set(gca,'XTick',0:30:360,'XTickLabel',[]); hold off;

    sp(3) = subplot(6,1,3);
    plot(rad2deg(t), acc*1000, 'k');        % acceleration
    hold on;
    plot(rad2deg(t(nt)), acc(nt)*1000,'ko','MarkerFaceColor','k');
    ylabel('acceleration [mm/s$^2$]'); axis tight;
    box on; grid on;
    set(gca,'XTick',0:30:360,'XTickLabel',[]); hold off;

    sp(4) = subplot(6,1,4);
    plot(rad2deg(t), jerk*1000, 'k');       % jerk
    hold on;
    plot(rad2deg(t(nt)), jerk(nt)*1000,'ko','MarkerFaceColor','k');
    ylabel('jerk [mm/s$^3$]'); axis tight;
    box on; grid on;
    set(gca,'XTick',0:30:360,'XTickLabel',[]); hold off;

    sp(5) = subplot(6,1,5);
    plot(rad2deg(t), rad2deg(phi), 'k');       % pressure angle
    hold on;
    plot(rad2deg(t(nt)), rad2deg(phi(nt)),'ko','MarkerFaceColor','k');

    tmp = [rad2deg(t) fliplr(rad2deg(t))];
    patch(tmp,[-30*ones(1,numel(t)) 30*ones(1,numel(t))],'b','FaceAlpha',0.2);
    ylabel('$\phi$ [degree]'); 
    xlabel('$\theta$ [degree]')
    axis tight; box on; grid on;
    set(gca,'XTick',0:30:360,'XTickLabel',0:30:360); hold off;

    sp(6) = subplot(6,1,6);

    plot(x_prime*1000,y_prime*1000,'b-'); hold on;      % prime circle
    plot(x_base*1000,y_base*1000,'b--');                % base circle
end

% Animate cam motion by rotating cam contour and ptich curve around cam center
rot_angle = t(nt);
R = [cos(rot_angle)     -sin(rot_angle); ...        % rotation matrix
     sin(rot_angle)     cos(rot_angle) ];

for n = 1:size(cam_in,2)
    cam_out(:,n) = R*cam_in(:,n);           % rotated cam contour
    pitch_out(:,n) = R*pitch_in(:,n);       % rotated pitch curve
end

patch(cam_out(1,:)*1000, cam_out(2,:)*1000,'r','FaceAlpha',0.5);
plot(pitch_out(1,:)*1000, pitch_out(2,:)*1000,'--k');

tt = deg2rad(0:360);
x_follower = r_follower*cos(tt);
y_follower = r_follower*sin(tt) + dis(nt) + r_prime;    % translate follower motion
patch(x_follower*1000,y_follower*1000,'k','FaceAlpha',0.5);

axis equal; box on; grid on;
xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('theta = %4.1f degree',rad2deg(t(nt))));

tmp = max(sqrt(x_pitch.^2 + y_pitch.^2)) + 2*r_follower;
axis([-1 1 -1 1]*tmp*1000);