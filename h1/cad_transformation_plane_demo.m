clear; close all; clc;

%% Load airplane stl file:
filename = 'aircraft.stl';
plane = stlread(filename);

figure(1); hold on;
trimesh(plane,'FaceColor','none','EdgeColor','k');
xlabel('x'); ylabel('y'); zlabel('z');
box on; grid on; view(3);
set(gca,'FontSize',14)

%% Assemble data into homogeneous coordinate:
points = [plane.Points'; ones(1,size(plane.Points,1))];

%% Define 3D transformation matricies in homogeneous coordinate
% Scaling matrix
S = @(sx, sy, sz) [sx, 0,  0, 0;
                    0, sy, 0, 0;
                    0, 0,  0, 1];

% Translation matrix
T = @(tx, ty, tz) [eye(4,3) [tx, ty, tz, 1]'];

% Rotation matrix around z axis: positive thetaz in degree counterclockwise
Rz = @(thetaz) [cosd(thetaz), -sind(thetaz), 0, 0;
                sind(thetaz),  cosd(thetaz), 0, 0;
                           0,             0, 1, 0;
                           0,             0, 0, 1];

% Rotation matrix around x axis: positive thetax in degree counterclockwise
Rx = @(thetax) [1,            0,             0, 0;
                0, cosd(thetax), -sind(thetax), 0;
                0, sind(thetax),  cosd(thetax), 0;
                0,            0,             0, 1];

% Rotation matrix around y axis: positive thetay in degree counterclockwise
Ry = @(thetay) [ cosd(thetay), 0, sind(thetay), 0;
                            0, 1,            0, 0;
                -sind(thetay), 0, cosd(thetay), 0;
                            0, 0,            0, 1];

%% Transform
transform_points = T(20,20,20)*Rx(45)*points;

transform_plane = triangulation(plane.ConnectivityList, ...
                 transform_points(1:3,:)'); % new points need to match Connectivity list

trimesh(transform_plane,'FaceColor','none','EdgeColor','r');

return;