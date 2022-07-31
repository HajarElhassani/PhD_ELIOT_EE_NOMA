clear all;
close all;
clc;

%% Parameter initialization
N=1e4; % number of parallel runs for the average regret
T=1e4; % number of iterations/training steps per run

% SNR thresholds (not rate thresholds)
threshold1_NOMA=1; % SNR threshold for user 1

Pmax=100; %power budget at the BS
Pc=1; %circuit power
var_h1=1; % channel variance of user 1
var_h2=1; % channel variance of user 2
sigma1=.1; % noise variance of user 1
sigma2=.1; % noise variance of user 2

% the different set of arms
arms0=[1 0.2;1 0.4;1 0.6;1 0.8;1 1;2 0.2;2 0.4;2 0.6;2 0.8;2 1];
arms1=[1 0.1;1 0.2;1 0.3;1 0.4;1 0.5;1 0.6;1 0.7;1 0.8;1 0.9;1 1;2 0.1;2 0.2;2 0.3;2 0.4;2 0.5;2 0.6;2 0.7;2 0.8;2 0.9;2 1];
arms2=[1 0.05;1 0.1;1 0.15;1 0.2;1 0.25;1 0.3;1 0.35;1 0.4;1 0.45;1 0.5;1 0.55;1 0.6;1 0.65;1 0.7;1 0.75;1 0.8;1 0.85;1 0.9;1 0.95;1 1;2 0.05;2 0.1;2 0.15;2 0.2;2 0.25;2 0.3;2 0.35;2 0.4;2 0.45;2 0.5;2 0.55;2 0.6;2 0.65;2 0.7;2 0.75;2 0.8;2 0.85;2 0.9;2 0.95;2 1];
% number of arms ion each arm set
a0=length(arms0);
a1=length(arms1);
a2=length(arms2);

%varying the SNR threshold of user 2
threshold2_NOMA=(1:.5:10);

% Parameters of learning algorithm UCB
alpha=.1; % upper bound parameter

%% compute the best fixed policy with CDIT for every QoS threshold 
for i=1:length(threshold2_NOMA)
    mu_best0= expected_offline_policy(arms0,Pmax,Pc,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    mu_best1= expected_offline_policy(arms1,Pmax,Pc,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    mu_best2= expected_offline_policy(arms2,Pmax,Pc,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    
    for n=1:N
        %% Initialization for UCB algorithm for every set of arms 
        number_of_sel_UCB0=ones(1,a0);% initial number of times arms were selected(the 6 arms were selected at least once)
        number_of_sel_UCB1=ones(1,a1); 
        number_of_sel_UCB2=ones(1,a2);
        
        empirical_mean_UCB0=zeros(1,a0);% initial reward sample (empirical reward) from each arm
        empirical_mean_UCB1=zeros(1,a1);
        empirical_mean_UCB2=zeros(1,a2);
        
        R_cum_UCB0=0;
        R_cum_UCB1=0;
        R_cum_UCB2=0;
        
        % initializing the iteration for each set of arms
        t=1;
        p=1;
        q=1;
        
        % initializing the regret for each set of arms  
        regret_UCB0=ones(1,t);
        regret_UCB1=ones(1,p);
        regret_UCB2=ones(1,q);
        
        %% compute the number of iterations required to reach 10% of the regret
        while(regret_UCB0(t)>0.1)
            %generate channel's gain
            G=channel_gain(1,sigma1,sigma2,var_h1,var_h2);               
            %UCB
            [cum_reward_UCB0(t),EE0,number_of_sel_UCB0,empirical_mean_UCB0,R_cum_UCB0]=UCB(t,arms0,alpha,Pmax,Pc,G,threshold1_NOMA,threshold2_NOMA(i),number_of_sel_UCB0,empirical_mean_UCB0,R_cum_UCB0);
            cum_regret_UCB0(t)=t*mu_best0-cum_reward_UCB0(t);
            regret_UCB0(t+1)=cum_regret_UCB0(t)/t;
            t=t+1;
        end
        
        while(regret_UCB1(p)>0.1)
            %generate channel's gain
            G=channel_gain(1,sigma1,sigma2,var_h1,var_h2);               
            %UCB
            [cum_reward_UCB1(p),EE1,number_of_sel_UCB1,empirical_mean_UCB1,R_cum_UCB1]=UCB(p,arms1,alpha,Pmax,Pc,G,threshold1_NOMA,threshold2_NOMA(i),number_of_sel_UCB1,empirical_mean_UCB1,R_cum_UCB1);
            cum_regret_UCB1(p)=p*mu_best1-cum_reward_UCB1(p);
            regret_UCB1(p+1)=cum_regret_UCB1(p)/p;
            p=p+1;
        end
        
        while(regret_UCB2(q)>0.1)
            %generate channel's gain
            G=channel_gain(1,sigma1,sigma2,var_h1,var_h2);               
            %UCB
            [cum_reward_UCB2(q),EE2,number_of_sel_UCB2,empirical_mean_UCB2,R_cum_UCB2]=UCB(q,arms2,alpha,Pmax,Pc,G,threshold1_NOMA,threshold2_NOMA(i),number_of_sel_UCB2,empirical_mean_UCB2,R_cum_UCB2);
            cum_regret_UCB2(q)=q*mu_best2-cum_reward_UCB2(q);
            regret_UCB2(q+1)=cum_regret_UCB2(q)/q;
            q=q+1;
        end
    
    % store the number of iterations for each run
    iteration0(n)=t-1;
    iteration1(n)=p-1;
    iteration2(n)=q-1; 
    
    end
    
    %averaging over all runs
    mean_iteration0(i)=abs(mean(iteration0));
    mean_iteration1(i)=abs(mean(iteration1));
    mean_iteration2(i)=abs(mean(iteration2));
   
end

%% plotting figures
figure(1)
semilogy(threshold2_NOMA,mean_iteration2,'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
semilogy(threshold2_NOMA,mean_iteration1,'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
semilogy(threshold2_NOMA,mean_iteration0,'-s','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
ylabel('Number of iterations');
xlabel('\Gamma_{2}^{th}');
legend('40 arms','20 arms','10 arms','location','best');
grid;