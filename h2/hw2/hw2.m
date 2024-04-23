clear all; close all; clc; format long;

name = 'Colleen Wang';
id = 'A16880246';
hw_num = 2;


%% Problem 1
% dimensions
a = 77; da = 0.5;
b = 30; db = 0.1;
c = 25; dc = 0.1;
g = 2; dg = 0.8;

% nominal value
d_nom = a - (b + c + g);

% worst case error
dd_wce = dg - da - db - dc;
p1a = [d_nom dd_wce];

% statistical error
sigma_g = dg/2; % desired
sigma_a = da/2; sigma_b = db/2; sigma_c = dc/2;
sigma_d = sqrt(sigma_g^2 - (sigma_a^2 + sigma_b^2 + sigma_c^2));
dd_se = sigma_d * 2;
p1b = [d_nom dd_se];


%% Problem 2
% parameters
c = 30;
db_lower = 0.3; db_upper = 0.5;
a = c - db_lower; % minimum
b = c + db_upper; % maximum

p2a = ( (a + b + c) / 3); % mean
p2b = c; % mode

x = 0.5;
cdf = @(x) (x <= c) .* ((x - a).^2 / ((b - a) * (c - a))) + ...
           (x > c) .* (1 - ((b - x).^2 / ((b - a) * (b - c))));
p2c = fzero(@(x) cdf(x) - 0.5, [a, b]);  % median

p2d = cdf(30); % probability of thickness < 30mm

cdf_mode = (c - a) / (b - a); % CDF at mode
inv_cdf = @(p) (p <= cdf_mode) .* (a + sqrt(p * (b - a) * (c - a))) + ...
               (p > cdf_mode) .* (b - sqrt((1 - p) * (b - a) * (b - c)));
p2e = inv_cdf(0.90); % 90th percentile


%% Problem 3
% parameters
S = 4; dS = 0.03;
p = 1.19; dp = 0.001; % rho
v = 33.4; dv = 0.5;
L = 135; dL = 0.5;

% nominal value
CL_nom = L / (0.5*p*v^2*S);
p3a = CL_nom; 

% wce
% CL_max = (L+dL) / ( 0.5 * (p-dp) * (v-dv)^2 * (S-dS) );
% CL_min = (L-dL) / ( 0.5 * (p+dp) * (v+dv)^2 * (S+dS) );
% p3b = max(abs(CL_nom - CL_max), abs(CL_nom - CL_min));
p3b = CL_nom*((dL/L) + 2*(dv/v) + (dp/p) + (dS/S)); 

% stat error
sigma_S = dS/2;
sigma_p = dp/2;
sigma_v = dv/2;
sigma_L = dL/2;
CL_se = CL_nom * sqrt( (sigma_L/L)^2 + (sigma_p/p)^2 + 2*(sigma_v/v)^2 + (sigma_S/S)^2 );
p3c = CL_se * 2;


% monte carlo analysis
rng('default'); % reset random number generator
n = 25000; % sample size
bins = 100;

L_samples = unifrnd(L-dL, L+dL, [n 1]); % uniform distribution
p_dist = makedist('Triangular', 'a', p-dp, 'b', p, 'c', p+dp);
p_samples = random(p_dist, [n 1]); % triangular distribution
v_dist = makedist('Triangular', 'a', v-dv, 'b', v, 'c', v+dv);
v_samples = random(v_dist, [n 1]);
S_samples = normrnd(S, dS / 2, [n 1]); % gaussian distribution for S

CL_samples = L_samples ./ (0.5 .* p_samples .* v_samples.^2 .* S_samples);

figure(1);
subplot(2, 3, 1);
histogram(L_samples, bins); grid on;
title('Histogram of Lift Force (L)');
xlabel('L (N)');
ylabel('Frequency');

subplot(2, 3, 2);
histogram(p_samples, bins); grid on;
title('Histogram of Fluid Density (ρ)');
xlabel('ρ (kg/m^3)');
ylabel('Frequency');

subplot(2, 3, 3);
histogram(v_samples, bins); grid on;
title('Histogram of Airspeed (v)');
xlabel('v (m/s)');
ylabel('Frequency');

subplot(2, 3, 4);
histogram(S_samples, bins); grid on;
title('Histogram of Airfoil Area (S)');
xlabel('S (m^2)');
ylabel('Frequency');

subplot(2, 3, 5);
histogram(CL_samples, bins); grid on;
title('Histogram of Lift Coefficient (C_L)');
xlabel('C_L');
ylabel('Frequency');

p3d = 'See figure 1';

% MLE
CL_gd = mle(CL_samples,'distribution','normal');
p3e = CL_gd(1); % mean
p3f = CL_gd(2); % standard deviation

p3g = normcdf(0.05, p3e, p3f);