clear all;
close all;
clc;

%% Parameter initialization
N=1e2; % number of runs
T=5000; % number of iterations/training steps



% SNR thresholds (not rate thresholds)
threshold1_NOMA=1; 
threshold2_NOMA=10;
% for OMA the SNR thresholds must be adjusted as the rate thershods must be
% identical
threshold1_OMA=threshold_conversion_OMA(threshold1_NOMA);
threshold2_OMA=threshold_conversion_OMA(threshold2_NOMA);

Pmax=100; %power budget
Pc=1; %circuit power
var_h1=1; %channel variance of link BS-user 1
var_h2=1; %channel variance of link BS-user 2
sigma1=.1;% noise variance of user 1
sigma2=.1;% noise variance of user 1

%set of arms
arms=[1 0.05;1 0.1;1 0.15;1 0.2;1 0.25;1 0.3;1 0.35;1 0.4;1 0.45;1 0.5;1 0.55;1 0.6;1 0.65;1 0.7;1 0.75;1 0.8;1 0.85;1 0.9;1 0.95;1 1;2 0.05;2 0.1;2 0.15;2 0.2;2 0.25;2 0.3;2 0.35;2 0.4;2 0.45;2 0.5;2 0.55;2 0.6;2 0.65;2 0.7;2 0.75;2 0.8;2 0.85;2 0.9;2 0.95;2 1];
a=length(arms); %number of arms
probability=1/a*ones(1,a); % uniform discrete distribution over the a arms


%UCB parameters
alpha=0.1; % upper bound parameter
%EXP3 parameters
gamma=min(1,sqrt(a*log(a)/((exp(1)-1)*T)));
eta=0.08;%gamma/a;


% offline fixed policy with CDIT      
mu_best= expected_offline_policy(arms,Pmax,Pc,threshold1_NOMA,threshold2_NOMA,sigma1,sigma2,var_h1,var_h2)
% offline fixed policy for OMA with CDIT  
mu_best_OMA= expected_offline_policy_OMA(arms,Pmax,Pc,threshold1_OMA,threshold2_OMA,sigma1,sigma2,var_h1,var_h2)

for n=1:N
    
    % generating channel gains for each run    
    G=channel_gain(T,sigma1,sigma2,var_h1,var_h2);
    
    %----------------------------------------------------------------------
    %% Initialization of learning algorithms
    %----------------------------------------------------------------------
    
    % Initialisation for UCB algorithm
    number_of_sel_UCB=ones(1,a); % initial number of times arms were selected(the 6 arms were selected at least once)
    empirical_mean_UCB=zeros(1,a);% initial reward sample (empirical reward) from each arm
    R_cum_UCB=0;
    
    % Initialisation for EXP3 algorithm
    number_of_sel_EXP3=ones(1,a);
    empirical_mean_EXP3=zeros(1,a);
    weights=ones(1,a);
    R_cum_EXP3=0;
 
    
    for t=1:T              
        %UCB
        [cum_reward_UCB(t),EE_UCB(t),number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB]=UCB(t,arms,alpha,Pmax,Pc,G(t,:),threshold1_NOMA,threshold2_NOMA,number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB);
                
        %EXP3        
        [cum_reward_EXP3(t),EE_EXP3(t),weights,number_of_sel_EXP3,empirical_mean_EXP3,R_cum_EXP3] = EXP3(t,arms,Pmax,Pc,weights,gamma,eta,G(t,:),threshold1_NOMA,threshold2_NOMA,number_of_sel_EXP3,empirical_mean_EXP3,R_cum_EXP3);          
    end
    
    %cumulative regret
    cum_regret_UCB(n,:)=(1:T)*mu_best-cum_reward_UCB;
    cum_regret_EXP3(n,:)=(1:T)*mu_best-cum_reward_EXP3;
    
    %regret
    regret_UCB(n,:)=cum_regret_UCB(n,:)./(1:T);
    regret_EXP3(n,:)=cum_regret_EXP3(n,:)./(1:T);
    
    %energy efficiency 
    EE_UCB_n(n,:)=EE_UCB;
    EE_EXP3_n(n,:)=EE_EXP3;
    
       
end

%averaging over all runs 
mean_cum_regret_UCB=mean(cum_regret_UCB);
mean_cum_regret_EXP3=mean(cum_regret_EXP3);

mean_regret_UCB=mean(regret_UCB);
mean_regret_EXP3=mean(regret_EXP3);

mean_EE_UCB=mean(EE_UCB_n);
mean_EE_EXP3=mean(EE_EXP3_n);


%% plotting figures

%cumulative regret
figure(1);
plot(1:T,mean_cum_regret_EXP3,'-<', 'MarkerIndices', 1:500:length(1:T));
hold on;
plot(1:T,mean_cum_regret_UCB,'->', 'MarkerIndices', 1:500:length(1:T));
ylabel('Cumulative regret');
xlabel('Iterations');
legend('EXP3','UCB','location','best');
grid on;

%regret
figure(2);
plot(1:T,mean_regret_EXP3,'-<', 'MarkerIndices', 1:500:length(1:T));
hold on;
plot(1:T,mean_regret_UCB,'->', 'MarkerIndices', 1:500:length(1:T));
ylabel('Regret');
xlabel('Iterations');
legend('EXP3','UCB','location','best');
grid on;

%energy efficiency
figure(3);
plot(1:T,mean_EE_EXP3,'-<','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 1:500:length(1:T));
hold on;
plot(1:T,(mu_best)*ones(1,T),'-*','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 1:500:length(1:T));
hold on;
plot(1:T,mean_EE_UCB,'->','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 1:500:length(1:T));
hold on;
plot(1:T,(mu_best_OMA)*ones(1,T),'--k','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 1:500:length(1:T));
ylabel('GEE (bits/J)');
xlabel('Iterations');
legend('EXP3','Best offline arm','UCB','OMA','location','best');
grid on;
