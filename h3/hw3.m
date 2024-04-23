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
    plot(x, y + offsets(i),'LineStyle','-','LineWidth',2);
end
xlabel('x'); ylabel('y');
title('Parabola with Offset Curves');
legend('y = x^2','D = 0.25','D = 0.5','D = 0.75','D=1','Location','northeast');
p1a = 'See figure 1';


dy = 2*x; dy2 = 2;
K = abs(dy2) ./ (1 + dy.^2).^(3/2);
R = 1 ./ K;

figure(2); hold on;
subplot(2, 1, 1);
plot(x, R, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('Radius of Curvature (R)');
legend('R');
title('Radius of y = x^2');

subplot(2, 1, 2);
plot(x, K, 'r-', 'LineWidth', 2);
xlabel('x');
ylabel('Curvature (K)');
legend('K');
title('Curvature of y = x^2');

p1b = 'See figure 2';

p1c = 'When the offset distance D exceeds the radius of curvature'