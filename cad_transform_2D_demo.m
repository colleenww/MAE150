clear; close all; clc;

%% Define parametric curve:
t = 0:360;
x = 16*sind(t).^3; % in degree
y = 13*cosd(t) - 5*cosd(2*t) - 2*cos(3*t) - cosd(4*t);

%% Assemble data into homogeneous coordinate:
points = [x; y; ones(1, length(x))];

%% Define an inline function for plotting
plot_points = @(list_of_points, cs) plot(list_of_points(1,:), list_of_points(2,:), [cs '-o']);

%% Define 2D transformation matricies in homogeneous coordinate
% scaling matrix
S = @(sx, sy) [sx, 0,  0; % anonymous function!
                0, sy, 0;
                0, 0,  1];

% Translation matrix
T = @(tx, ty) [eye(3,2) [tx ty 1]']; % creates diagonal matrix

% Rotation matrix: positive theta in degree counterclockwise
R = @(theta) [cosd(theta), -sind(theta), 0;
              sind(theta),  cosd(theta), 0;
              0,            0,          1];

%% Transform data

figure(1); hold on;
plot_points(points, 'k');
axis equal;
xlabel('x'); ylabel('y');
box on; grid on;
set(gca,'FontSize',14)

transform_points = T(30,30)*points;
plot_points(transform_points,'r');