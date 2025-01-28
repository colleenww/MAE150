clear all;   % This clears all workspaces
close all;   % This closes all figures 
clc;         % This clears the command window
format long; 

%% Student name, ID number, and homework_number are required
name = 'Colleen Wang';
id = 'A16880246';
hw_num = 5;

%% Problem 1
E = 200e9;              % Young's Modulus (Pa)
A = pi*(0.05/2)^2;      % cross-sectional area (m^2)
L = 2;
n_node = 7;             % number of nodes
n_element = 11;         % number of elements
DOF = 2;                % degree of freedom per node
scale_factor = 1000;    % used to amplify nodal displacement in figure
font_size = 15;         % font size used in figure

% create mesh
n = 1;
node(n).pos = [0,0];    % x-, y- coordinates
node(n).dis = [0,0];    % x-, y- components
node(n).force = [NaN, NaN];  % x-, y- components: NaN = unknown

n = 2;
node(n).pos = [L,0];
node(n).dis = [NaN, NaN];
node(n).force = [0, -15e3];

n = 3;
node(n).pos = [2*L,0];
node(n).dis = [NaN, NaN];
node(n).force = [0, -15e3];

n = 4;
node(n).pos = [3*L,0];
node(n).dis = [NaN, 0];
node(n).force = [0, NaN];

n = 5;
node(n).pos = [L/2, 1];
node(n).dis = [NaN, NaN];
node(n).force = [0, 0];

n = 6;
node(n).pos = [3*L/2, 1];
node(n).dis = [NaN, NaN];
node(n).force = [0, 0];

n = 7;
node(n).pos = [5*L/2, 1];
node(n).dis = [NaN, NaN];
node(n).force = [15e3*cosd(30), 15e3*sind(30)];

% Elements
n = 1;
element(n).node = [1,2]; % smaller node number comes first
element(n).E = E;
element(n).A = A;

n = 2;
element(n).node = [2,3]; 
element(n).E = E;
element(n).A = A;

n = 3;
element(n).node = [3,4];
element(n).E = E;
element(n).A = A;

n = 4;
element(n).node = [1,5];
element(n).E = E;
element(n).A = A;

n = 5;
element(n).node = [2,5];
element(n).E = E;
element(n).A = A;

n = 6;
element(n).node = [2,6];
element(n).E = E;
element(n).A = A;

n = 7;
element(n).node = [3,6];
element(n).E = E;
element(n).A = A;

n = 8;
element(n).node = [3,7];
element(n).E = E;
element(n).A = A;

n = 9;
element(n).node = [4,7];
element(n).E = E;
element(n).A = A;

n = 10;
element(n).node = [5,6];
element(n).E = E;
element(n).A = A;

n = 11;
element(n).node = [6,7];
element(n).E = E;
element(n).A = A;

p1a = zeros(4, 4, n_element);

% compute elemental stiffness matrix
K_e_fun = @(E,A,L,beta) E*A/L *...
    [cosd(beta)^2,          cosd(beta)*sind(beta),      -cosd(beta)^2,          -cosd(beta)*sind(beta);
     cosd(beta)*sind(beta), sind(beta)^2,               -cosd(beta)*sind(beta), -sind(beta)^2,
     -cosd(beta)^2,         -cosd(beta)*sind(beta),     cosd(beta)^2,           cosd(beta)*sind(beta);
     -cosd(beta)*sind(beta),-sind(beta)^2,              cosd(beta)*sind(beta),  sind(beta)^2];

for n = 1:numel(element)
    % compute length L
    dx = (node(element(n).node(2)).pos(1) - node(element(n).node(1)).pos(1));
    dy = (node(element(n).node(2)).pos(2) - node(element(n).node(1)).pos(2));
    element(n).L = sqrt(dx^2 + dy^2);

    % obtain slope angle beta
    beta = atand(dy/dx);
    if beta < 0; beta = beta + 180; end
    if dy < 0; beta = beta + 180; end
    element(n).beta = beta;

    % elemental stiffness matrix
    K_e = K_e_fun(element(n).E, element(n).A,element(n).L,element(n).beta);
    element(n).K_e = K_e;
    p1a(:,:,n) = K_e;

    % expand elemental stiffness to global dimension
    K_global_e = zeros(DOF * numel(node));
    iL = element(n).node(1);
    iH = element(n).node(2);
    K_global_e([2*iL-1, 2*iL, 2*iH-1, 2*iH],[2*iL-1, 2*iL, 2*iH-1, 2*iH]) = K_e;
    element(n).K_global_e = K_global_e;
end

% Assemble global stiffness matrix
K_global = zeros(DOF*numel(node));
for n = 1:numel(element)
    K_global = K_global + element(n).K_global_e;
end
p1b = K_global;

% Apply boundary conditions to reduce the equilibrium equations
% identify row and column to keep: where the displacement is unknown (i.e. equal to NaN)
U_global = [node.dis]';
row_col_to_keep = find(isnan(U_global));

% reduce global stiffness matrix
K_reduced = K_global(row_col_to_keep, row_col_to_keep);
p1c = K_reduced;

% reduce force vector
F_global = [node.force]';
F_reduced = F_global(row_col_to_keep);

% solve equilibrium equations
% solve for reduced global displacement vector
U_reduced = K_reduced \ F_reduced;

% update global displacement vector
U_global(row_col_to_keep) = U_reduced;

% update global force vector
F_global = K_global*U_global;

% pass updated global displacement and force back to nodes
for n = 1:numel(node)
    node(n).dis = U_global(2*n-1:2*n)';
    node(n).force = F_global(2*n-1:2*n)';
end
p1d = reshape(U_global, DOF, n_node)';
p1e = reshape(F_global, DOF, n_node)';

% plot original and deformed trusses
figure('unit','inches','Position',[0.5 0.5 0 5]);
hold on;

% element
for n = 1:numel(element)
    x_elem = [node(element(n).node(1)).pos(1), node(element(n).node(2)).pos(1)];
    y_elem = [node(element(n).node(1)).pos(2), node(element(n).node(2)).pos(2)];
    plot(x_elem,y_elem,'--ko','MarkerFaceColor','k','LineWidth',3);

    dx_elem = [node(element(n).node(1)).dis(1), node(element(n).node(2)).dis(1)];
    dy_elem = [node(element(n).node(1)).dis(2), node(element(n).node(2)).dis(2)];
    p1 = patch(x_elem + scale_factor * dx_elem, ...
               y_elem + scale_factor * dy_elem, ...
               sqrt(dx_elem.^2 + dy_elem.^2),'EdgeColor','interp','LineWidth',3);
    text(mean(x_elem),mean(y_elem), sprintf('[%d]',n), 'Color','b','FontSize',font_size);
end
colorbar; colormap('jet');

% node and force vectors
for n = 1:numel(node)
    xd = node(n).pos(1);
    yd = node(n).pos(2);
    max_force = max(abs(F_global));
    q = quiver(xd,yd,node(n).force(1)/max_force, node(n).force(2)/max_force,0.2,'Color','r','LineWidth',2,'MaxHeadSize',2);
    text(xd, yd, sprintf('%d',n),'Color','b','FontSize',font_size);
end
box on; grid on; axis equal;
xlabel('x [m]');
ylabel('y [m]');
title(sprintf('Truss deformation using scaling factor of %.0f',scale_factor));
set(gca,'FontSize',font_size);
p1f = 'See figure 1';


%% Problem 2
r = 0.02;  % Radius in meters
I = (pi * r^4) / 4;  % Moment of inertia for a circular cross section
L = 1;  % Length of the beam in meters
P = 1000;  % Applied load in N

% nodes and elements
nodes = [0, 0.1, 0.2, 0.8, 0.9, L];
num_nodes = length(nodes);
element_lengths = diff(nodes);
num_elements = length(element_lengths);

% 3d stiffness matrix
p2a = zeros(4, 4, num_elements);
for i = 1:num_elements
    Le = element_lengths(i);
    k = (E * I) / (Le^3) * ...
        [12,  6*Le, -12,  6*Le;
         6*Le, 4*Le^2, -6*Le, 2*Le^2;
        -12, -6*Le,  12, -6*Le;
         6*Le, 2*Le^2, -6*Le, 4*Le^2];
    p2a(:,:,i) = k;
end

% global stiffness matrix
K_global = zeros(2*num_nodes, 2*num_nodes);
for i = 1:num_elements
    K_e = p2a(:,:,i);
    node_indices = 2*(i-1) + (1:4); % Mapping local DOFs to global DOFs
    K_global(node_indices, node_indices) = K_global(node_indices, node_indices) + K_e;
end
p2b = K_global;

% reduced global stiffness matrix
K_reduced = K_global(3:end, 3:end); % beam clamped at x = 0
p2c = K_reduced;

% load vector
F = zeros(2*num_nodes, 1);
F(end-1) = P;
F_reduced = F(3:end);

U_reduced = K_reduced \ F_reduced;

% full displacement vector
U = zeros(2*num_nodes, 1);
U(3:end) = U_reduced;

% deflections (w) and bending angles (theta)
w = U(1:2:end);  % displacements
theta = U(2:2:end);  % rotations
p2d = [w, theta];

% force (F) and moment (M)
F_M_global = K_global * U;
F_nodes = F_M_global(1:2:end);
M_nodes = F_M_global(2:2:end);
p2e = [F_nodes, M_nodes];

% Theoretical deflection calculation
w_theory = @(x) -P .* x.^2 / (6 * E * I) .* (3 * L - x);  % negative for downward deflection
dx = 0.1; % grid resolution
x_theory = 0:dx:L;
w_theory_values = w_theory(x_theory);

w_interpolated = interp1(nodes, -w, x_theory, 'pchip');  % 'pchip' for piecewise cubic Hermite interpolating polynomial

figure(2); hold on;
plot(x_theory, w_theory_values, 'r-', 'DisplayName', 'Theoretical Solution');
plot(x_theory, w_interpolated, 'b--', 'DisplayName', 'FEM Analysis');
plot(nodes, -w, 'ko', 'MarkerSize', 6, 'DisplayName', 'Nodes');
xlabel('Position Along the Beam (m)');
ylabel('Deflection (m)');
title('Comparison of FEM vs. Theoretical Deflection Analysis');
legend('show');
grid on;
p2f = 'See figure 2';


%% Problem 4
d_i = 0.05; % Initial diameter in meters (5 cm)
L = 2; % Length of elements in meters
P1_i = 10e3; % Initial load P1 in Newtons (10 kN)
P2_i = 10e3; % Initial load P2 in Newtons (10 kN)
A_i = pi * (d_i / 2)^2; % Initial cross-sectional area

% Doubling the loads
P1_new = 2 * P1_i;
P2_new = 2 * P2_i;

k_initial = E * A_i / L;
delta = (P1_i + P2_i) * L / (A_i * E);
A_new = (P1_new + P2_new) * L / (delta * E);
d_new = 2 * sqrt(A_new / pi);
d_new_mm = round(d_new * 1000); % convert to mm and round
p4a = d_new_mm*10^-3;
p4a = 0.07358;