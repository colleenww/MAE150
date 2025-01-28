clear all;   % This clears all workspaces
close all;   % This closes all figures 
clc;         % This clears the command window
format long; 

%% Student name, ID number, and homework_number are required
name = 'Colleen Wang';
id = 'A16880246';
hw_num = 6;

%% Problem 1
E = 200e9;
s = 0.3;
I = s^4 / 12;
L_beam = 12;
n_node = 4;
n_element = n_node - 1;
DOF = 2;
font_size = 12;
plot_theory = true;         % flag to include theoretical solution

%% Symbolically define shape function
syms x L N1 N2 N3 N4
N1 = 1 - 3*(x/L)^2 + 2*(x/L)^3;
N2 = x - 2*x^2/L + x^3/L^2;
N3 = 3*x^2/L^2 - 2*x^3/L^3;
N4 = -x^2/L + x^3/L^2;

%% Symbolically define load and theoretical solution
syms q w_theory theta_theory moment_theory shear_theory

q = -100;      % distributed load
P = -1;        % point load
w_theory = P*x^2/6/E/I * (3*L_beam - x);
theta_theory = P*x/2/E/I * (2*L_beam - x);
moment_theory = P*(L_beam - x);
shear_theory = -P;

%% Create mesh
% Nodes
n = 1;
node(n).pos = 0;            % x-coordinates along the beam
node(n).dis = [0, 0];       % [deflection, angle]
node(n).force_p = [0, 0];   % [point force, point moment]

n = 2;
node(n).pos = 4;
node(n).dis = [NaN, NaN];
node(n).force_p = [-1000, 0];

n = 3;
node(n).pos = 8;
node(n).dis = [NaN, NaN];
node(n).force_p = [500, 0];

n = 4;
node(n).pos = 12;
node(n).dis = [NaN, NaN];
node(n).force_p = [-2000, 0];

% Elements
for n = 1:n_element
    element(n).node = [n, n+1];
    element(n).E = E;
    element(n).I = I;
end

%% Convert distributed load to equivalent nodal forces
for n = 1:numel(element)
    eL = abs(node(element(n).node(2)).pos(1) - node(element(n).node(1)).pos(1));
    if n <= 2 % Apply distributed load only to the first two elements
        element(n).force_q1 = [ eval(subs(int(q*N1,x,[0,eL]),eL)), ...      % F1_q
                                eval(subs(int(q*N2,x,[0,eL]),eL)) ];        % M1_q
        element(n).force_q2 = [ eval(subs(int(q*N3,x,[0,eL]),eL)), ...      % F2_q
                                eval(subs(int(q*N4,x,[0,eL]),eL)) ];        % M2_q
    else
        element(n).force_q1 = [0, 0];
        element(n).force_q2 = [0, 0];
    end
end

%% Combine point and distributed load at each node
for n = 1:numel(node)
    node(n).force_pq = node(n).force_p;
    for e = 1:numel(element)
        if element(e).node(1) == n
            node(n).force_pq = node(n).force_pq + element(e).force_q1;
        elseif element(e).node(2) == n
            node(n).force_pq = node(n).force_pq + element(e).force_q2;
        end
    end
end

%% Compute elemental stiffness matrix
% Beam elemental stiffness matrix
K_e_fun = @(E,I,L) E*I/L^3 * ...
                            [12,  6*L,  -12,  6*L; ...
                             6*L, 4*L^2, -6*L, 2*L^2; ...
                            -12, -6*L,   12,  -6*L; ...
                             6*L, 2*L^2, -6*L, 4*L^2];

for n = 1:numel(element)
    % Compute element length L
    element(n).eL = abs( node(element(n).node(2)).pos(1) ...
                       - node(element(n).node(1)).pos(1));

    % Elemental stiffness matrix:
    K_e = K_e_fun(element(n).E, element(n).I, element(n).eL);
    element(n).K_e = K_e;

    % Expand elemental stiffness to global dimension
    K_global_e = zeros(DOF * numel(node));
    iL = element(n).node(1);
    iH = element(n).node(2);
    K_global_e([2*iL-1, 2*iL, 2*iH-1, 2*iH], [2*iL-1, 2*iL, 2*iH-1, 2*iH]) = K_e;
    element(n).K_global_e = K_global_e;
end

%% Assemble global stiffness matrix
K_global = zeros(DOF*numel(node));
for n = 1:numel(element)
    K_global = K_global + element(n).K_global_e;
end
p1a = K_global;

%% Apply boundary conditions to reduce the equilibrium equations
U_global = [node.dis]';
row_col_to_keep = find(isnan(U_global));

% Reduce global stiffness matrix
K_reduced = K_global(row_col_to_keep, row_col_to_keep);
p1b = K_reduced;

% Reduce force vector
F_global = [node.force_pq]';
F_reduced = F_global(row_col_to_keep);

%% Solve equilibrium equations
% Solve for reduced global displacement vector
U_reduced = K_reduced \ F_reduced;

% Update global displacement vector
U_global(row_col_to_keep) = U_reduced;

% Update global force vector
F_global = K_global * U_global;

%% Pass updated global displacement and force back to nodes
% Also compute the reaction force
for n = 1:numel(node)
    node(n).dis = U_global(2*n-1:2*n)';
    node(n).force_global = F_global(2*n-1:2*n)';
    node(n).force_reaction = node(n).force_global - node(n).force_pq;
end

% Extract reaction force and moment at the wall (node 1)
reaction_force = node(1).force_reaction(1);
reaction_moment = node(1).force_reaction(2);

% Store the reaction force and moment in a row vector
p1d = [reaction_force, reaction_moment];

%% Store displacement vectors in a 4x2 matrix
p1c = reshape([node.dis], 2, n_node)';

%% Compute w(x), theta(x), and V(x) for each element
for n = 1:numel(element)
    x_local = linspace(node(element(n).node(1)).pos(1), ...
                       node(element(n).node(2)).pos(1), 100);
    element(n).x_local = x_local;

    w_local = N1*node(element(n).node(1)).dis(1) ...
            + N2*node(element(n).node(1)).dis(2) ...
            + N3*node(element(n).node(2)).dis(1) ...
            + N4*node(element(n).node(2)).dis(2);
    element(n).w_local = eval(subs(w_local,{x,L}, ...
                    {x_local-x_local(1),x_local(end)-x_local(1)}));

    theta_local = diff(N1,x)*node(element(n).node(1)).dis(1) ...
                + diff(N2,x)*node(element(n).node(1)).dis(2) ...
                + diff(N3,x)*node(element(n).node(2)).dis(1) ...
                + diff(N4,x)*node(element(n).node(2)).dis(2);
    element(n).theta_local = eval(subs(theta_local,{x,L}, ...
                    {x_local-x_local(1),x_local(end)-x_local(1)}));

    shear_local_expr = E * I * ...
                 ( diff(N1,x,3)*node(element(n).node(1)).dis(1) ...
                 + diff(N2,x,3)*node(element(n).node(1)).dis(2) ...
                 + diff(N3,x,3)*node(element(n).node(2)).dis(1) ...
                 + diff(N4,x,3)*node(element(n).node(2)).dis(2) );
    element(n).shear_local = eval(subs(shear_local_expr,{x,L}, ...
                {x_local-x_local(1),x_local(end)-x_local(1)}));

    moment_local_expr = E * I * ...
                  ( diff(N1,x,2)*node(element(n).node(1)).dis(1) ...
                  + diff(N2,x,2)*node(element(n).node(1)).dis(2) ...
                  + diff(N3,x,2)*node(element(n).node(2)).dis(1) ...
                  + diff(N4,x,2)*node(element(n).node(2)).dis(2) );
    element(n).moment_local = eval(subs(moment_local_expr,{x,L}, ...
                    {x_local-x_local(1),abs(x_local(end)-x_local(1))}));
end

%% Plot results
figure('unit','in','Position',[0.5 0.5 0 0]);

% Applied loads
% sp(1) = subplot(5,1,1); hold on;
% plot([0 L_beam],[0 0],'b-');        % reference zero line
% for n = 1:numel(node)               % point force
%     quiver(node(n).pos,0,0,node(n).force_p(1),0.1,...
%         'r','filled','LineWidth',3);
% end
% 
% for n = 1:numel(node)
%     if node(n).force_p(2) > 0       % positive point moment: green
%        plot(node(n).pos,0,'go',...,
%             'MarkerFaceColor','g','MarkerSize',10);
%     elseif node(n).force_p(2) < 0   % negative point moment: magenta
%         plot(node(n).pos,0,'mo',...,
%             'MarkerFaceColor','m','MarkerSize',10);
%     end
% end
% 
% for n = 1:numel(element)-1            % distributed force
%     bar(element(n).x_local(1:10:end),...
%         eval(subs(q,x,element(n).x_local(1:10:end))),'b','FaceAlpha',0.5);
% end
% ylabel({'applied load','not-to-scale'});
% xlim([-0.05 1.05]* L_beam); box on; grid on;
% set(gca,'YTickLabel',[],'FontSize',font_size);

% Deflection
sp(1) = subplot(4,1,1); hold on;
for n = 1:numel(element)
    plot(element(n).x_local, element(n).w_local,'LineWidth',8);
end
% if plot_theory
%     xl = linspace(0,L_beam,100);
%     plot(xl,eval(subs(w_theory,x,xl)),'--k');
% end
ylabel('deflection [m]'); xlim([-0.05 1.05]* L_beam); box on; grid on;
set(gca,'FontSize',font_size);

% Bending angle theta(x)
sp(2) = subplot(4,1,2); hold on;
for n = 1:numel(element)
    plot(element(n).x_local, element(n).theta_local,'LineWidth',8);
end
% if plot_theory
%     xl = linspace(0,L_beam,100);
%     plot(xl,eval(subs(theta_theory,x,xl)),'--k');
% end
ylabel('slope'); xlim([-0.05 1.05]* L_beam); box on; grid on;
set(gca,'FontSize',font_size);

% Bending moment diagram N(x)
sp(3) = subplot(4,1,3); hold on;
for n = 1:numel(element)
    plot(element(n).x_local, element(n).moment_local,'LineWidth',8);
end
% if plot_theory
%     plot(xl,eval(subs(moment_theory,x,xl)),'--k');
% end
ylabel('moment [N m]'); xlim([-0.05 1.05]* L_beam); box on; grid on;
set(gca,'FontSize',font_size);

% Shear force diagram V(x)
sp(4) = subplot(4,1,4); hold on;
for n = 1:numel(element)
    plot(element(n).x_local, element(n).shear_local,'LineWidth',8);
end
% if plot_theory
%     plot(xl,eval(subs(shear_theory,x,xl)),'--k');
% end
xlabel('distance along the beam x [m]');
ylabel('shear [N]'); xlim([-0.05 1.05]* L_beam); box on; grid on;
sgtitle('FEA of beam');
set(gca,'FontSize',font_size);
p1f = 'See figure 1';

%% Problem 2
ansys_dvalues = [0, -1.37e-3, -4.6e-3, -8.76e-3];
node_positions = [0, 4, 8, 12];

figure(2); hold on;
for n = 1:numel(element)
    plot(element(n).x_local, element(n).w_local,'LineWidth',2);
end
plot(node_positions, ansys_dvalues, '-ko','LineWidth',2);
xlabel('distance along the beam x [m]')
ylabel('deflection [m]'); xlim([-0.05 1.05]* L_beam); box on; grid on;
title('Comparison of Ansys vs. MATLAB Deflection Values');
legend('MATLAB','Ansys');
p2b = 'See figure 2';

%% Problem 3
rng('default');

% Parameters
n = 500; % Number of samples
nom_s = 0.3; % Nominal side length in meters (30 cm)
tol_s = 0.01; % Tolerance in meters (1 cm)
E = 200e9; % Elastic modulus in Pascals
L_beam = 12; % Beam length in meters
P2 = -1000; % Load at 4 m
P3 = 500; % Load at 8 m
P4 = -2000; % Load at 12 m
q = -100; % Distributed load

% Generate samples for s from a normal distribution
s_samples = normrnd(nom_s, tol_s, [1, n]);

% Pre-allocate array for deflection results
deflections = zeros(1, n);

% Loop through each sample
for i = 1:n
    % Update moment of inertia for current sample
    s_current = s_samples(i);
    I_current = s_current^4 / 12;

    % Node positions and initial displacements and forces
    node(1).pos = 0;  node(1).dis = [0, 0];  node(1).force_p = [0, 0];
    node(2).pos = 4;  node(2).dis = [NaN, NaN];  node(2).force_p = [P2, 0];
    node(3).pos = 8;  node(3).dis = [NaN, NaN];  node(3).force_p = [P3, 0];
    node(4).pos = 12; node(4).dis = [NaN, NaN];  node(4).force_p = [P4, 0];

    % Element definitions
    for n = 1:3
        element(n).node = [n, n+1];
        element(n).E = E;
        element(n).I = I_current;
        element(n).eL = abs(node(element(n).node(2)).pos - node(element(n).node(1)).pos);
    end

    % Shape functions
    syms x L N1 N2 N3 N4
    N1 = 1 - 3*(x/L)^2 + 2*(x/L)^3;
    N2 = x - 2*x^2/L + x^3/L^2;
    N3 = 3*x^2/L^2 - 2*x^3/L^3;
    N4 = -x^2/L + x^3/L^2;

    % Convert distributed load to equivalent nodal forces
    for n = 1:numel(element)
        eL = element(n).eL;
        if n <= 2
            element(n).force_q1 = [ eval(subs(int(q*N1,x,[0,eL]),eL)), eval(subs(int(q*N2,x,[0,eL]),eL)) ];
            element(n).force_q2 = [ eval(subs(int(q*N3,x,[0,eL]),eL)), eval(subs(int(q*N4,x,[0,eL]),eL)) ];
        else
            element(n).force_q1 = [0, 0];
            element(n).force_q2 = [0, 0];
        end
    end

    % Combine point and distributed load at each node
    for n = 1:numel(node)
        node(n).force_pq = node(n).force_p;
        for e = 1:numel(element)
            if element(e).node(1) == n
                node(n).force_pq = node(n).force_pq + element(e).force_q1;
            elseif element(e).node(2) == n
                node(n).force_pq = node(n).force_pq + element(e).force_q2;
            end
        end
    end

    % Elemental stiffness matrix function
    K_e_fun = @(E,I,L) E*I/L^3 * [12,  6*L,  -12,  6*L; 6*L, 4*L^2, -6*L, 2*L^2; -12, -6*L, 12, -6*L; 6*L, 2*L^2, -6*L, 4*L^2];

    % Assemble global stiffness matrix
    K_global = zeros(8, 8);
    for n = 1:numel(element)
        K_e = K_e_fun(element(n).E, element(n).I, element(n).eL);
        iL = element(n).node(1);
        iH = element(n).node(2);
        indices = [2*iL-1, 2*iL, 2*iH-1, 2*iH];
        K_global(indices, indices) = K_global(indices, indices) + K_e;
    end

    % Apply boundary conditions
    U_global = [node.dis]';
    row_col_to_keep = find(isnan(U_global));

    % Reduce global stiffness matrix and force vector
    K_reduced = K_global(row_col_to_keep, row_col_to_keep);
    F_global = [node.force_pq]';
    F_reduced = F_global(row_col_to_keep);

    % Solve for reduced global displacement vector
    U_reduced = K_reduced \ F_reduced;

    % Update global displacement vector
    U_global(row_col_to_keep) = U_reduced;

    % Extract deflection at the free end (node 4)
    deflections(i) = U_global(7); % deflection at node 4
end

% Plot histograms of the side of the cross-section s and the deflection w
figure;
subplot(1, 2, 1);
histogram(s_samples, 10, 'Normalization', 'pdf');
title('Histogram of Side of Cross Section s');
xlabel('s (m)');
ylabel('Probability Density');
grid on;

subplot(1, 2, 2);
histogram(deflections, 10, 'Normalization', 'pdf');
title('Histogram of Deflection w');
xlabel('Deflection (m)');
ylabel('Probability Density');
grid on;

p3a = 'See figure 3';

% Fit a Gaussian distribution to the deflection data using MLE
pd_deflection = fitdist(deflections', 'Normal');
deflection_mean = pd_deflection.mu;
deflection_std_dev = pd_deflection.sigma;

p3b = deflection_mean;
p3c = deflection_std_dev;

% Compute tolerance for the deflection based on 95% confidence interval
confidence_level = 0.95;
z_score = norminv(0.5 + confidence_level / 2);
deflection_tolerance = z_score * deflection_std_dev;

p3d = deflection_tolerance;