clear; close all; clc;
rng('default'); % reset random seed

n = 1e5;

%% Setup bucket part A:
nom_A = 30;
tol_A = 2;
distribution_type = 'uniform';
lower = nom_A - tol_A;
upper = nom_A + tol_A;
pd_A = makedist(distribution_type,'lower',lower,'upper',upper);
bucket_A = random(pd_A, [1,n]);

%% Setup bucket part B
nom_B = 40;
tol_B = 3;
distribution_type = 'uniform';
lower = nom_B - tol_B;
upper = nom_B + tol_B;
pd_B = makedist(distribution_type,'lower',lower,'upper',upper);
bucket_B = random(pd_B,[1,n]);

%% Setup bucket part C
nom_C = 70;
tol_C = 2;
distribution_type = 'uniform';
lower = nom_C - tol_C;
upper = nom_C + tol_C;
pd_C = makedist(distribution_type,'lower',lower,'upper',upper);
bucket_C = random(pd_C,[1,n]);

%% Compute bucket gap;
bucket_gap = bucket_C - bucket_A - bucket_B;

%% Fit normal distribution to bucket gap using MLE method
pd_gap = mle(bucket_gap,'distribution','normal'); % 2 elements: mu and sigma
mu_gap = pd_gap(1);
sigma_gap = pd_gap(2);
mle_fit = @(x) normpdf(x,mu_gap,sigma_gap); % distribution of the bucket_Gap

P_gap_it_0mm = normcdf(0,mu_gap,sigma_gap)
P_gap_gt_2mm = 1-normcdf(2,mu_gap,sigma_gap)

%% Plot results:
figure('Units','inches','Position',[1 1 15 10]);

sp(1) = subplot(3,2,1); hold on;
histogram(bucket_A,'BinWidth',0.1,'Normalization','pdf');
plot(mean(bucket_A)*[1,1],[0,1],'--r','LineWidth',1);
plot(median(bucket_A)*[1,1],[0,1],'--c','LineWidth',1);
text(min(bucket_A),0.95,sprintf('Mean: %5.3f',mean(bucket_A)),'FontSize',14,'Color','r');
text(min(bucket_A),0.9,sprintf('Mode: %5.3f',my_mode(bucket_A)),'FontSize',14,'Color','b');
text(min(bucket_A),0.85,sprintf('Median: %5.3f',median(bucket_A)),'FontSize',14,'Color','c');
text(min(bucket_A),0.8,sprintf('St. dev.: %5.3f',std(bucket_A)),'FontSize',14);
ylabel('PDF'); xlabel('Part A size (mm)');
box on; grid on;
set(gca,'FontSize',12);
