clear all;
close all;
clc;

%% Parameter initialization
N=1e2; % number of parallel runs for the average regret
T=1e4; % number of iterations/training steps per run

% SNR thresholds (not rate thresholds)
threshold1_NOMA=1; % SNR threshold for user 1

Pmax=100; %power budget at the BS
var_h1=1; % channel variance of user 1
var_h2=0.1; % channel variance of user 2
sigma1=1; % noise variance of user 1
sigma2=1; % noise variance of user 2


% the different set of arms
arms0=[1 0.25;2 0.25];
arms1=[1 0.125;1 0.25;1 0.375;2 0.125;2 0.25;2 0.375];
arms2=[1 0.0625;1 0.125;1 0.1875;1 0.25;1 0.3125;1 0.375;1 0.4375;2 0.0625;2 0.125;2 0.1875;2 0.25;2 0.3125;2 0.375;2 0.4375];
arms3=[1 0.03125;1 0.0625;1 0.09375;1 0.125;1 0.15625;1 0.1875;1 0.21875;1 0.25;1 0.28125;1 0.3125;1 0.34375;1 0.375;1 0.40625;1 0.4375;1 0.46875;2 0.03125;2 0.0625;2 0.09375;2 0.125;2 0.15625;2 0.1875;2 0.21875;2 0.25;2 0.28125;2 0.3125;2 0.34375;2 0.375;2 0.40625;2 0.4375;2 0.46875];
% number of arms ion each arm set
a0=length(arms0);
a1=length(arms1);
a2=length(arms2);
a3=length(arms3);

%varying the SNR threshold of user 2
threshold2_NOMA=(1:20);

% Parameters of learning algorithm UCB
alpha=1; % upper bound parameter


%% compute the best fixed policy with CDIT for every QoS threshold 
for i=1:length(threshold2_NOMA)
    [mu_best0,mu_worst0]= expected_offline_policy(arms0,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    [mu_best1,mu_worst1]= expected_offline_policy(arms1,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    [mu_best2,mu_worst2]= expected_offline_policy(arms2,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    [mu_best3,mu_worst3]= expected_offline_policy(arms3,Pmax,threshold1_NOMA,threshold2_NOMA(i),sigma1,sigma2,var_h1,var_h2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% iteration for N runs 
    for n=1:N
        
        %% Initialization for UCB algorithm for every set of arms 
        
        number_of_sel_UCB0=ones(1,a0);% initial number of times arms were selected
        number_of_sel_UCB1=ones(1,a1); 
        number_of_sel_UCB2=ones(1,a2);
        number_of_sel_UCB3=ones(1,a3);
        
        empirical_mean_UCB0=zeros(1,a0);% initial reward sample (empirical reward) from each arm
        empirical_mean_UCB1=zeros(1,a1);
        empirical_mean_UCB2=zeros(1,a2);
        empirical_mean_UCB3=zeros(1,a3);
        
        R_cum_UCB0=0;
        R_cum_UCB1=0;
        R_cum_UCB2=0;
        R_cum_UCB3=0;
        
        % initializing the iteration for each set of arms
        t=1;
        p=1;
        q=1;
        r=1;
        
        % initializing the regret for each set of arms      
        regret_UCB0=ones(1,t);
        regret_UCB1=ones(1,p);
        regret_UCB2=ones(1,q);
        regret_UCB3=ones(1,r);
        
        %% compute the number of iterations required to reach 10% of the regret
        while(regret_UCB0(t)>0.1)
            %generate channels gain
            G=channel_gain(1,sigma1,sigma2,var_h1,var_h2);               
            % UCB algorithm
            [cum_reward_UCB0(t),outage_UCB0,number_of_sel_UCB0,empirical_mean_UCB0,R_cum_UCB0]=UCB(t,arms0,alpha,Pmax,G,threshold1_NOMA,threshold2_NOMA(i),number_of_sel_UCB0,empirical_mean_UCB0,R_cum_UCB0);
            cum_regret_UCB0(t)=t*mu_best0-cum_reward_UCB0(t);
            % update the regret
            regret_UCB0(t+1)=cum_regret_UCB0(t)/t;
            t=t+1;
        end
        
        while(regret_UCB1(p)>0.1)
            %generate channels gain
            G=channel_gain(1,sigma1,sigma2,var_h1,var_h2);               
            % UCB algorithm
            [cum_reward_UCB1(p),outage_UCB1,number_of_sel_UCB1,empirical_mean_UCB1,R_cum_UCB1]=UCB(p,arms1,alpha,Pmax,G,threshold1_NOMA,threshold2_NOMA(i),number_of_sel_UCB1,empirical_mean_UCB1,R_cum_UCB1);
            cum_regret_UCB1(p)=p*mu_best1-cum_reward_UCB1(p);
            % update the regret
            regret_UCB1(p+1)=cum_regret_UCB1(p)/p;
            p=p+1;
        end
        
        while(regret_UCB2(q)>0.1)
            %generate channels gain
            G=channel_gain(1,sigma1,sigma2,var_h1,var_h2);               
            % UCB algorithm
            [cum_reward_UCB2(q),outage_UCB2,number_of_sel_UCB2,empirical_mean_UCB2,R_cum_UCB2]=UCB(q,arms2,alpha,Pmax,G,threshold1_NOMA,threshold2_NOMA(i),number_of_sel_UCB2,empirical_mean_UCB2,R_cum_UCB2);
            cum_regret_UCB2(q)=q*mu_best2-cum_reward_UCB2(q);
            % update the regret
            regret_UCB2(q+1)=cum_regret_UCB2(q)/q;
            q=q+1;
        end
        
        while(regret_UCB3(r)>0.1)
            %generate channels gain
            G=channel_gain(1,sigma1,sigma2,var_h1,var_h2);               
            % UCB algorithm
            [cum_reward_UCB3(r),outage_UCB3,number_of_sel_UCB3,empirical_mean_UCB3,R_cum_UCB3]=UCB(r,arms3,alpha,Pmax,G,threshold1_NOMA,threshold2_NOMA(i),number_of_sel_UCB3,empirical_mean_UCB3,R_cum_UCB3);
            cum_regret_UCB3(r)=r*mu_best3-cum_reward_UCB3(r);
            % update the regret
            regret_UCB3(r+1)=cum_regret_UCB3(r)/r;
            r=r+1;
        end
    
    % store the number of iterations for each run
    iteration0(n)=t-1;
    iteration1(n)=p-1;
    iteration2(n)=q-1; 
    iteration3(n)=r-1;
    
    end
    
    % averaging over all runs
    mean_iteration0(i)=abs(mean(iteration0));
    mean_iteration1(i)=abs(mean(iteration1));
    mean_iteration2(i)=abs(mean(iteration2));
    mean_iteration3(i)=abs(mean(iteration3));
   
end


%% plotting figures
figure(1)
semilogy(threshold2_NOMA,mean_iteration3,'-o','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
semilogy(threshold2_NOMA,mean_iteration2,'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
semilogy(threshold2_NOMA,mean_iteration1,'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices',1:2:length(threshold2_NOMA));
hold on;
semilogy(threshold2_NOMA,mean_iteration0,'-s','MarkerSize',7,'LineWidth',1.3);
ylabel('Number of iterations');
xlabel('\Gamma_{2}^{th}');
legend('30 arms','14 arms','6 arms','2 arms','location','best');
grid;