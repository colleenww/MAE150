%% MAE 150 Homework Solution Tempate
% 
% Please name your file hwX.m where the X is the homework number.
%  For example, homework 1 should be named hw1.m
%
% All homework assignments need to be uploaded to CANVAS by midnight on
%   on Friday the week they are due.  
%
% Make sure that your homework script runs and all answers are specified 
%   before submission.
%
% If your script does not run, or your answers are not reported in the 
%  correct format, you will be given no credit.

%% The following commands are required at the very top of the homework file
clear all;   % This clears all workspaces
close all;   % This closes all figures 
clc;         % This clears the command window
format long; 

%% Student name, ID number, and homework_number are required
name = 'Colleen Wang';
id = 'A16880246';
hw_num = 1;

    
%% Programming problems: 
%  Include all code in your script. Answers are to be reported as pX= where 
%  X is the problem number.
%  For example the answer to problem1 should be reported as p1=
%  if it is a single part question, and p1a= if it is a multipart
%  question. 

% This is how to report your answer for a single part question
% p1 = 100;
% 
% % This is how to report your answer for a multi-part question
% p2a = 'test';
% p2b = 2234;
% p2c = 11/5;
% p2d = floor(p3c);
% p2e = [1 4 5 7 1e2];

%% Problem 1
data = readtable('US_COVID19.txt', 'Delimiter', '\t', 'Headerlines', 9);

p1a = data.PostiveResults(data.Date == 20200704);
p1b = data.PostiveResults(data.Date == 20201126);

figure(1); hold on;
bar(data.YearDay, data.PostiveResults, 'b');
plot(data.YearDay(data.Date == 20200704), p1a, 'cd', 'MarkerFaceColor', 'c', 'MarkerSize', 10);
plot(data.YearDay(data.Date == 20201126), p1b, 'mo', 'MarkerFaceColor', 'm', 'MarkerSize', 10);
running_avg = movmean(data.PostiveResults, [7 7]);
running_avg(1:7) = NaN;
running_avg(end-6:end) = NaN;
plot(data.YearDay, running_avg, 'r-', 'LineWidth', 1);
xlabel('Year Day'); ylabel('Number of Positive Cases'); legend('Positive Cases', 'July 4th', 'Thanksgiving', '15-day Running Average', 'Location', 'northeast');
title('Number of Positive COVID-19 Cases in the US');
p1c = 'See figure 1';

%% Problem 3
% Define parametric curve:
t = 0:360;
x = cosd(t);
y = zeros(1,length(x));
z = sind(t).*cosd(t);

% Assemble data into homogeneous coordinate:
points = [x; y; z; ones(1, length(x)) ];

% Define an inline function for plotting
plot_points = @(list_of_points, cs) ...
    plot3(list_of_points(1,:), list_of_points(2,:), ...
    list_of_points(3,:), [cs '-']);

% Rotation matrix around z axis: positive thetaz in degree counterclockwise
Rz = @(thetaz) [cosd(thetaz), -sind(thetaz), 0, 0;
                sind(thetaz),  cosd(thetaz), 0, 0;
                           0,             0, 1, 0;
                           0,             0, 0, 1];

% Transform data
figure(2); hold on;
plot_points(points, '');
xlabel('x'); ylabel('y'); zlabel('z'); 
title('Curve Rotated About the z-axis');
box on; grid on;
view(3); axis equal;
set(gca,'FontSize',14);

for n = 5:5:175
    transform_points = Rz(n)*points;
    plot_points(transform_points,'');
end 

p3 = 'See figure 2';

%% Problem 4
filename = 'aircraft.stl';
plane = stlread(filename);

figure(3); hold on;
trimesh(plane,'FaceColor','none','EdgeColor','k'); hold on;
xlabel('x'); ylabel('y'); zlabel('z');
box on; grid on; view(3); axis equal;
set(gca,'FontSize',14)
title('Aircraft Movement in 10s');

% Assemble data into homogeneous coordinate:
points = [plane.Points'; ones(1,size(plane.Points,1))];

% Define 3D transformation matricies in homogeneous coordinate
% Translation matrix
T = @(tx, ty, tz) [eye(4,3) [tx, ty, tz, 1]'];

% Rotation matrix around x axis: positive thetax in degree counterclockwise
Rx = @(thetax) [1,            0,             0, 0;
                0, cos(thetax), -sin(thetax), 0;
                0, sin(thetax),  cos(thetax), 0;
                0,            0,             0, 1];

% movement parameters
speed_x = -2; % m/s
angular_velocity_x = 0.1; % counterclockwise
time = 10; % seconds

displacement_x = speed_x * time;
rotation_x = angular_velocity_x * time;

% Transform
transform_points = T(displacement_x,0,0)*Rx(rotation_x)*points;

transform_plane = triangulation(plane.ConnectivityList, ...
                 transform_points(1:3,:)'); % new points need to match Connectivity list

trimesh(transform_plane,'FaceColor','none','EdgeColor','r');
legend('Initial Position','Final Position','Location','southeast');

highest_point_z = max(transform_points(3,:));
highest_point_index = find(transform_points(3,:) == highest_point_z);
p4a = transform_points(1:3, highest_point_index);
p4b = 'See figure 3';

return;