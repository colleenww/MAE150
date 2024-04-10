clear all; close all; clc;

%% Define parametric curve:
t = 0:360;
x = 16*sind(t).^3;
y = zeros(1,length(x));
z = 13*cosd(t) - 5*cosd(2*t) - 2*cosd(3*t) - cosd(4*t);

%% Assemble data into homogeneous coordinate:
points = [x; y; z; ones(1, length(x)) ];

%% Define an inline function for plotting
plot_points = @(list_of_points, cs) ...
    plot3(list_of_points(1,:), list_of_points(2,:), ...
    list_of_points(3,:), [cs '-o']);

%% Define 3d transformation matricies in homogeneous coordinate
% Scaling matrix
S = @(sx, sy, sz) [sx, 0,  0, 0;
                    0, sy, 0, 0;
                    0, 0, sz, 0;
                    0, 0,  0, 1];

% Translation matrix
T = @(tx, ty, tz) [eye(4,3) [tx, ty, tz, 1]'];

% Rotation matrix around z axis: positive thetaz in degree counterclockwise
Rz = @(thetaz) [cosd(thetaz), -sind(thetaz), 0, 0;
                sind(thetaz),  cosd(thetaz), 0, 0;
                           0,             0, 1, 0;
                           0,             0, 0, 1];

%% Transform data
figure(1); hold on;
plot_points(points,'k');
xlabel('x'); ylabel('y'); zlabel('z');
box on; grid on;
view(3); axis([-1 1 -1 1 -1 1]*10,'equal');
set(gca,'FontSize',14)

% % to make sure the code works, perform a simple transformation
% transform_points = T(30,0,30)*points;
% plot_points(transform_points,'r');

% % extrusion
% for n = 1:10
%     % translate the heart in the y direction 10 times
%     transform_points = T(0,n,0)*points;
%     plot_points(transform_points,'r');
% end

% revolution
for n = 10:10:350
    transform_points = Rz(n)*points;
    plot_points(transform_points,'r');
end 


