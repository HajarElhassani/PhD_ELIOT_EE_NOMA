clear all;
close all;
clc;
 
%% parameters initialization
Pmax=100; % power budget at the BS
var_h1=1; % channel variance at user 1
var_h2=0.1; % channel variance at user 2
sigma1=1; % noise variance at user 1
sigma2=1; % noise variance at user 2
% for OMA the SNR thresholds must be adjusted as the rate thershods must be
% identical

% SNR thresholds (not rate thresholds)
threshold1_NOMA=1;
threshold1_OMA=threshold_conversion_OMA(threshold1_NOMA);

% set of arms 
arms1=[1 0.25;2 0.25];
arms2=[1 0.125;1 0.25;1 0.375;2 0.125;2 0.25;2 0.375];
arms3=[1 0.0625;1 0.125;1 0.1875;1 0.25;1 0.3125;1 0.375;1 0.4375;2 0.0625;2 0.125;2 0.1875;2 0.25;2 0.3125;2 0.375;2 0.4375];
arms4=[1 0.03125;1 0.0625;1 0.09375;1 0.125;1 0.15625;1 0.1875;1 0.21875;1 0.25;1 0.28125;1 0.3125;1 0.34375;1 0.375;1 0.40625;1 0.4375;1 0.46875;2 0.03125;2 0.0625;2 0.09375;2 0.125;2 0.15625;2 0.1875;2 0.21875;2 0.25;2 0.28125;2 0.3125;2 0.34375;2 0.375;2 0.40625;2 0.4375;2 0.46875];
% who performs SIC
users=[1 2];

%offline fixed policy with different QoS thresholds
threshold2_NOMA=(1:20);
alpha=(0:1e-5:0.5);

%% compute the best offline fixed policy with CDIT for each threshold and different sets of arms

for i=1:length(threshold2_NOMA)
    % for OMA the SNR thresholds must be adjusted as the rate thershods must be identical
    threshold2_OMA(i)=threshold_conversion_OMA(threshold2_NOMA(i));   
    
    %% outage of OMA:
    mu_OMA=expected_value_OMA(threshold1_OMA,threshold2_OMA(i),Pmax,sigma1,sigma2,var_h1,var_h2);
    outage_OMA(i)=1-mu_OMA;
    
    % find the best fixed policy (alpha: the fraction of power) with CDIT by exhaustive search 
    for j=1:length(alpha)
        [mu_best,mu_worst]= expected_offline_policy_alpha(users,alpha(j),Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
        outage_NOMA=1-mu_best;
        %compute the outage gap between OMA and NOMA
        gap(j)=100*(outage_OMA(i)-outage_NOMA)/outage_OMA(i); 
    end
    % take the alpha (fraction of power) which gives the maximum gap between OMA and NOMA
    [mxv_max,idx_max]=max(gap);
    alpha_max = alpha(idx_max);
    [mxv_min,idx_min]=min(gap);
    alpha_min = alpha(idx_min);
   
    %% outage of NOMA
    % the best fixed policy with CDIT by exhaustive search 
    [mu_best_max,mu_worst]= expected_offline_policy_alpha(users,alpha_max,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    outage_NOMA_max(i)=1-mu_best_max;

    
    % the best fixed policy with CDIT for the set arms 1
    [mu_best1,mu_worst1]= expected_offline_policy(arms1,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    outage_NOMA1(i)=1-mu_best1;
   
    % the best fixed policy with CDIT for the set arms 2
    [mu_best2,mu_worst2]= expected_offline_policy(arms2,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    outage_NOMA2(i)=1-mu_best2;
    
    % the best fixed policy with CDIT for the set arms 3
    [mu_best3,mu_worst3]= expected_offline_policy(arms3,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    outage_NOMA3(i)=1-mu_best3;
    
    % the best fixed policy with CDIT for the set arms 4
    [mu_best4,mu_worst4]= expected_offline_policy(arms4,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    outage_NOMA4(i)=1-mu_best4;
    
end

%% plot outage of OMA and outages of different arms
figure(1)
plot(threshold2_NOMA,100*outage_NOMA4,'-o','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,100*outage_NOMA3,'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,100*outage_NOMA2,'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,100*outage_NOMA1,'-s','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,100*outage_NOMA_max,'->','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
plot(threshold2_NOMA,100*outage_OMA,'--k','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
xlabel('\Gamma_2^{th}');
ylabel('Outage (%)');
grid on;
legend('30 arms','14 arms','6 arms','2 arms','Optimal solution','OMA','location','best');

