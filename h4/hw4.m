clear all;   % This clears all workspaces
close all;   % This closes all figures 
clc;         % This clears the command window
format long; 

%% Student name, ID number, and homework_number are required
name = 'Colleen Wang';
id = 'A16880246';
hw_num = 4;

%% Problem 2
k1 = 20; k2 = 40; k3 = 120; f2 = 10; f3 = -20;
k_12 = k1 + k2;
K = [k_12, -k_12, 0; -k_12, k_12+k3, -k3; 0, -k3, k3];
p2a = K;
K_reduced = K(2:end, 2:end);
F = [0; f2; f3];
F_reduced = [f2; f3];
d_reduced = K_reduced \ F_reduced;
d = [0; d_reduced];
p2b = d;
F = K * d;
p2c = F;

%% Problem 3
E = 210e9;
d = 0.05;
A = pi * (d/2)^2;

node1 = [1, 0];
node2 = [0, 1];

L = sqrt((node2(1) - node1(1))^2 + (node2(2) - node1(2))^2);
cos_theta = (node2(1) - node1(1)) / L;
sin_theta = (node2(2) - node1(2)) / L;

% Local stiffness matrix (2x2)
k_local = (E * A / L) * [1, -1; -1, 1];

% Transformation matrix (4x4)
T = [cos_theta, sin_theta, 0, 0; 
     -sin_theta, cos_theta, 0, 0;
     0, 0, cos_theta, sin_theta; 
     0, 0, -sin_theta, cos_theta];

% Global stiffness matrix (4x4)
k_global = T' * blkdiag(k_local, k_local) * T;

% Display or return k_global
disp(k_global);


% Correct transformation matrix (4x4)
T = [cos_theta, sin_theta, 0, 0; 
    -sin_theta, cos_theta, 0, 0;
    0, 0, cos_theta, sin_theta; 
    0, 0, -sin_theta, cos_theta];

% Global stiffness matrix (4x4) using proper 4x4 block diagonal
k_global = T' * blkdiag(k_local, k_local) * T;

% Display or return k_global
disp(k_global);
