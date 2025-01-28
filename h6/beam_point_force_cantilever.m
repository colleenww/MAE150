clear all;   % This clears all workspaces
close all;   % This closes all figures 
clc;         % This clears the command window
format long; 

%% Set up parameters
E = 200e9;
I = 1.186e-4;
L_beam = 7.5;
n_node = 2;
n_element = n_node - 1;
DOF = 2;
font_size = 12;
plot_theory = true;         % flag to include theoretical solution

%% Symbolically define shpae function
syms x L N1 N2 N3 N4        % convert q(x)
N1 = 1 - 3*(x/L)^2 + 2*(x/L)^3;
N2 = x - 2*x^2/L + x^3/L^2;
N3 = 3*x^2/L^2 - 2*x^3/L^3;
N4 = -x^2/L + x^3/L^2;

%% Symbolically define load and theoretical solution
syms q w_theory theta_theory moment_theory shear_theory

% cantilever with point load at the other end
q = 0;      % distributed load
P = -1;     % point load
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
node(n).pos = 7.5;
node(n).dis = [NaN, NaN];
node(n).force_p = [-1, 0];

% Elements
for n = 1:n_element
    element(n).node = [n, n+1];
    element(n).E = E;
    element(n).I = I;
end

%% Convert distributed load to point load for each element
for n = 1:numel(element)
    eL = abs(node(element(n).node(2)).pos(1) ...
        - node(element(n).node(1)).pos(1));
    element(n).force_q1 = [ eval(subs(int(q*N1,x,[0,eL]),eL)), ...      % F1_q
                            eval(subs(int(q*N2,x,[0,eL]),eL)) ];        % M1_q
    element(n).force_q2 = [ eval(subs(int(q*N3,x,[0,eL]),eL)), ...      % F2_q
                            eval(subs(int(q*N4,x,[0,eL]),eL)) ];        % M2_q
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
                       - node(element(n).node(1)).pos(1) );

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

%% Apply boundary conditions to reduce the equilibrim equations
U_global = [node.dis]';
row_col_to_keep = find(isnan(U_global));

% Reduce global stiffness matrix
K_reduced = K_global(row_col_to_keep, row_col_to_keep);

% Reduce force vector
F_global = [node.force_pq]';
F_reduced = F_global(row_col_to_keep);

%% Solve equilibrium equations
% Solve for reduced global displacement vector
U_reduced = K_reduced \ F_reduced;

% Update global displacement vector
U_global(row_col_to_keep) = U_reduced;

% Update global force vector
F_global = K_global*U_global;

%% Pass updated global displacement and force back to nodes
% Also compute the reaction force
for n = 1:numel(node)
    node(n).dis = U_global(2*n-1:2*n)';
    node(n).force_global = F_global(2*n-1:2*n)';
    node(n).force_reaction = node(n).force_global - node(n).force_pq;
end

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
sp(1) = subplot(5,1,1); hold on;
plot([0 L_beam],[0 0],'b-');        % reference zero line
for n = 1:numel(node)               % point force
    quiver(node(n).pos,0,0,node(n).force_p(1),0.1,...
        'r','filled','LineWidth',3);
end

for n = 1:numel(node)
    if node(n).force_p(2) > 0       % positive point moment: green
       plot(node(n).pos,0,'go',...,
            'MarkerFaceColor','g','MarkerSize',10);
    elseif node(n).force_p(2) < 0   % negative point moment: magenta
        plot(node(n).pos,0,'mo',...,
            'MarkerFaceColor','m','MarkerSize',10);
    end
end

for n = 1:numel(element)            % distributed force
    bar(element(n).x_local(1:10:end),...
        eval(subs(q,x,element(n).x_local(1:10:end))),'b','FaceAlpha',0.5);
end
ylabel({'applied load','not-to-scale'});
xlim([-0.05 1.05]* L_beam); box on; grid on;
set(gca,'YTickLabel',[],'FontSize',font_size);

% Deflection
sp(2) = subplot(5,1,2); hold on;
for n = 1:numel(element)
    plot(element(n).x_local, element(n).w_local,'LineWidth',8);
end
if plot_theory
    xl = linspace(0,L_beam,100);
    plot(xl,eval(subs(w_theory,x,xl)),'--k');
end
ylabel('deflection [m]'); xlim([-0.05 1.05]* L_beam); box on; grid on;
set(gca,'FontSize',font_size);

% Bending angle theta(x)
sp(3) = subplot(5,1,3); hold on;
for n = 1*numel(element)
    plot(element(n).x_local, element(n).theta_local,'LineWidth',8);
end
if plot_theory
    plot(xl,eval(subs(theta_theory,x,xl)),'--k');
end
ylabel('slope'); xlim([-0.05 1.05]* L_beam); box on; grid on;
set(gca,'FontSize',font_size);

% Bending moment diagram N(x)
sp(4) = subplot(5,1,4); hold on;
for n = 1:numel(element)
    plot(element(n).x_local, element(n).moment_local,'LineWidth',8);
end
if plot_theory
    plot(xl,eval(subs(moment_theory,x,xl)),'--k');
end
ylabel('moment [N m]'); xlim([-0.05 1.05]* L_beam); box on; grid on;
set(gca,'FontSize',font_size);

% Shear force diagram V(x)
sp(5) = subplot(5,1,5); hold on;
for n = 1:numel(element)
    plot(element(n).x_local, element(n).shear_local,'LineWidth',8);
end
if plot_theory
    plot(xl,eval(subs(shear_theory,x,xl)),'--k');
end
xlabel('distance along the beam x [m]');
ylabel('shear [N]'); xlim([-0.05 1.05]* L_beam); box on; grid on;
sgtitle('FEA of beam');
set(gca,'FontSize',font_size);
