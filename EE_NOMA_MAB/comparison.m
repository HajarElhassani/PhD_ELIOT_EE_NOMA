clear all;
close all;
clc;
 
%% parameters initialization
Pmax=100; %power budget
Pc=1; %circuit power
var_h1=1;% channel variance at user 1
var_h2=1;% channel variance at user 2
sigma1=.1;% noise variance at user 1
sigma2=.1;% noise variance at user 2
% for OMA the SNR thresholds must be adjusted as the rate thershods must be
% identical

% SNR thresholds (not rate thresholds)
threshold1_NOMA=1;
threshold1_OMA=threshold_conversion_OMA(threshold1_NOMA);

% set of arms 
arms1=[1 0.2;1 0.4;1 0.6;1 0.8;1 1;2 0.2;2 0.4;2 0.6;2 0.8;2 1];
arms2=[1 0.1;1 0.2;1 0.3;1 0.4;1 0.5;1 0.6;1 0.7;1 0.8;1 0.9;1 1;2 0.1;2 0.2;2 0.3;2 0.4;2 0.5;2 0.6;2 0.7;2 0.8;2 0.9;2 1];
arms3=[1 0.05;1 0.1;1 0.15;1 0.2;1 0.25;1 0.3;1 0.35;1 0.4;1 0.45;1 0.5;1 0.55;1 0.6;1 0.65;1 0.7;1 0.75;1 0.8;1 0.85;1 0.9;1 0.95;1 1;2 0.05;2 0.1;2 0.15;2 0.2;2 0.25;2 0.3;2 0.35;2 0.4;2 0.45;2 0.5;2 0.55;2 0.6;2 0.65;2 0.7;2 0.75;2 0.8;2 0.85;2 0.9;2 0.95;2 1];
% who performs SIC
users=[1 2];

%offline fixed policy with different QoS thresholds
threshold2_NOMA=(1:10);
beta=(1e-2:1e-2:1);

%% compute the best offline fixed policy with CDIT for each threshold and different sets of arms

for i=1:length(threshold2_NOMA)
    % for OMA the SNR thresholds must be adjusted as the rate thershods must be identical   
    threshold2_OMA(i)=threshold_conversion_OMA(threshold2_NOMA(i));   
    
    %% find the suboptimal fixed policy with CDIT by exhaustive search on beta
    %first compute all reward values for each beta
    for j=1:length(beta)
        EE_NOMA(j)= offline_policy_beta(users,beta(j),Pmax,Pc,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    end
    % take beta with best value of the reward 
    [mxv_max,idx_max]=max(EE_NOMA);
    beta_max = beta(idx_max);
    % get the suboptimal fixed policy with CDIT corresponding to the best taken beta
    EE_NOMA_max(i)= EE_NOMA(idx_max)
       
    % get the optimal fixed policy by exhaustive search
    EE_NOMA_op(i) = optimal_NOMA(threshold1_NOMA,threshold2_NOMA(i),Pmax,Pc,sigma1,sigma2,var_h1,var_h2);
       
    % get the suboptimal fixed policy with CDIT for arms 1
    EE_NOMA1(i)= expected_offline_policy(arms1,Pmax,Pc,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
   
    % get the suboptimal fixed policy with CDIT for arms 2
    EE_NOMA2(i)= expected_offline_policy(arms2,Pmax,Pc,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    
    % get the suboptimal fixed policy with CDIT for arms 3
    EE_NOMA3(i)= expected_offline_policy(arms3,Pmax,Pc,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    
    
end

%% plot 
figure(1)
plot(threshold2_NOMA,EE_NOMA3,'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,EE_NOMA2,'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,EE_NOMA1,'-s','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,EE_NOMA_max,'->','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,EE_NOMA_op,'-<','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:length(threshold2_NOMA));
xlabel('\Gamma_2^{th}');
ylabel('GEE (bits/J)');
grid on;
legend('40arms','20 arms','10 arms','suboptimal','optimal','location','best');

