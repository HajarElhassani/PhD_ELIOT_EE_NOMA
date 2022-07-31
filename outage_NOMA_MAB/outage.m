clear all;
close all;
clc;
%% Parameter initialization
N=1e2; % number of parallel runs for the average regret
T=1e4; % number of iterations/training steps per run


% SNR thresholds (not rate thresholds)
threshold1_NOMA=1; 
threshold2_NOMA=3;
% for OMA the SNR thresholds must be adjusted as the rate thershods must be
% identical
threshold1_OMA=threshold_conversion_OMA(threshold1_NOMA);
threshold2_OMA=threshold_conversion_OMA(threshold2_NOMA);

Pmax=100; % power budget at the base station

var_h1=1; % channel variance of link BS-user 1
var_h2=0.1; % channel variance of link BS-user 2

sigma1=1; % noise variance of user 1
sigma2=1; % noise variance of user 2


% set of arms 
arms = [1 0.0625;1 0.125;1 0.1875;1 0.25;1 0.3125;1 0.375;1 0.4375;2 0.0625;2 0.125;2 0.1875;2 0.25;2 0.3125;2 0.375;2 0.4375];

%number of arms
a=length(arms);

% uniform discrete distribution over the a arms
probability=1/a*ones(1,a);


% EXP3 learning parameters
gamma=min(1,sqrt(a*log(a)/((exp(1)-1)*T)));
eta=0.02; % gamma/a;

% Parameters of learning algorithms 
alpha=0.5; % upper bound parameter


% offline fixed policy with CDIT    
[mu_best,mu_worst]= expected_offline_policy(arms,Pmax,threshold1_NOMA,threshold2_NOMA,sigma1,sigma2,var_h1,var_h2);

% offline OMA with CDIT
mu_OMA = expected_value_OMA(threshold1_OMA,threshold2_OMA,Pmax,sigma1,sigma2,var_h1,var_h2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iteration for N runs 
for n=1:N
    
    % generating channel gains for each run 
    G=channel_gain(T,sigma1,sigma2,var_h1,var_h2);

    %----------------------------------------------------------------------
    %% Initialization of learning algorithms
    %----------------------------------------------------------------------
   
    
    % Initialisation for epsilon greedy algorithm
    number_of_sel_epsilon=zeros(1,a);
    empirical_mean_epsilon=zeros(1,a);
    R_cum_epsilon=0;
    
    % Initialisation for EXP3 algorithm
    number_of_sel_EXP3=zeros(1,a);
    weights=ones(1,a);
    R_cum_EXP3=0;
    
        
    % Initialization for UCB algorithm
    number_of_sel_UCB=zeros(1,a); % initial number of times arms were selected
    empirical_mean_UCB=zeros(1,a);% initial rewards sample (empirical reward) from each arm
    R_cum_UCB=0;  
    

    %----------------------------------------------------------------------
    %% EXP3 
    %----------------------------------------------------------------------
    for t=1:T
             
        % call one iteration of EXP3        
        [cum_reward_EXP3(t), outage_EXP3(t), weights,number_of_sel_EXP3,R_cum_EXP3] = EXP3(t,gamma,eta,arms,Pmax,weights,G(t,:),threshold1_NOMA,threshold2_NOMA,number_of_sel_EXP3,R_cum_EXP3);
      
    end
    %% epsilon-greedy, UCB
    % Initialization
    % test each arm of index t for t = 1 : a
    %----------------------------------------------------------------------
    for t = 1 : a
        
        
        % call one iteration of epsilon-greedy init 
        %[cum_reward_epsilon(t),outage_epsilon(t),number_of_sel_epsilon,empirical_mean_epsilon,R_cum_epsilon] = epsilon_greedy_init(t,arms,Pmax,G(t,:),threshold1_NOMA,threshold2_NOMA,vepsilon,probability,empirical_mean_epsilon,number_of_sel_epsilon,R_cum_epsilon);
       
        % call one iteration of UCB init
        [cum_reward_UCB(t),outage_UCB(t),number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB]=UCB_init(t,arms,Pmax,G(t,:),threshold1_NOMA,threshold2_NOMA,number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB);
    
    
    
    end
    %----------------------------------------------------------------------
    %% epsilon-greedy, UCB
    % generic iterations for t > a
    %----------------------------------------------------------------------
    for t = a+1 : T
        
        
        % call one generic iteration of epsilon-greedy
        %[cum_reward_epsilon(t),outage_epsilon(t),number_of_sel_epsilon,empirical_mean_epsilon,R_cum_epsilon] = epsilon_greedy(t,arms,Pmax,G(t,:),threshold1_NOMA,threshold2_NOMA,vepsilon,probability,empirical_mean_epsilon,number_of_sel_epsilon,R_cum_epsilon);
       
        
        % call one generic iteration of UCB
        [cum_reward_UCB(t),outage_UCB(t),number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB]=UCB(t,arms,alpha,Pmax,G(t,:),threshold1_NOMA,threshold2_NOMA,number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %cumulative regret
    cum_regret_UCB(n,:)=(1:T)*mu_best-cum_reward_UCB;
    cum_regret_EXP3(n,:)=(1:T)*mu_best-cum_reward_EXP3;
    %cum_regret_epsilon(n,:)=(1:T)*mu_best-cum_reward_epsilon;

    
    %regret
    regret_UCB(n,:)=cum_regret_UCB(n,:)./(1:T);
    regret_EXP3(n,:)=cum_regret_EXP3(n,:)./(1:T);
    %regret_epsilon(n,:)=cum_regret_epsilon(n,:)./(1:T);

    %outage 
    outage_UCB_n(n,:)=outage_UCB;
    outage_EXP3_n(n,:)=outage_EXP3;
    %outage_epsilon_n(n,:)=outage_epsilon;

       
end

%% averaging over all runs 
mean_cum_regret_UCB=mean(cum_regret_UCB);
mean_cum_regret_EXP3=mean(cum_regret_EXP3);
%mean_cum_regret_epsilon=mean(cum_regret_epsilon);


mean_regret_UCB=mean(regret_UCB);
mean_regret_EXP3=mean(regret_EXP3);
%mean_regret_epsilon=mean(regret_epsilon);

mean_outage_UCB=mean(outage_UCB_n);
mean_outage_EXP3=mean(outage_EXP3_n);
%mean_outage_epsilon=mean(outage_epsilon_n);


%--------------------------------------------------------------------------
%% Plotting figures
%--------------------------------------------------------------------------

%cumulative regret
figure(1);
%plot(1:T,mean_cum_regret_epsilon,'-d', 'LineWidth',1.3,'MarkerIndices', 100:T/10:T);
%hold on;
plot(1:T,mean_cum_regret_EXP3,'-<', 'LineWidth',1.3,'MarkerIndices', 100:T/10:T);
hold on;
plot(1:T,mean_cum_regret_UCB,'->', 'LineWidth',1.3,'MarkerIndices', 100:T/10:T);
ylabel('Cumulative regret');
xlabel('Iterations');
legend('EXP3','UCB');
grid on;

%regret
figure(2);
%plot(1:T,mean_regret_epsilon,'-d','LineWidth',1.3, 'MarkerIndices', 100:T/10:T);
%hold on;
plot(1:T,mean_regret_EXP3,'-<', 'LineWidth',1.3,'MarkerIndices', 100:T/10:T);
hold on;
plot(1:T,mean_regret_UCB,'->', 'LineWidth',1.3,'MarkerIndices', 100:T/10:T);
ylabel('Regret');
xlabel('Iterations');
legend('EXP3','UCB');
grid on;

%outage
figure(3);
%plot(1:T,100*mean_outage_epsilon,'-d', 'LineWidth',1.3,'MarkerIndices', 100:T/10:T);
%hold on
plot(1:T,100*mean_outage_EXP3,'-<', 'MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 100:T/10:T);
hold on
plot(1:T,100*(1-mu_best)*ones(1,T),'-*','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 100:T/10:T);
hold on
plot(1:T,100*mean_outage_UCB,'->','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 100:T/10:T);
hold on
plot(1:T,100*(1-mu_OMA)*ones(1,T),'--k','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 100:T/10:T);
ylabel('Outage (%)');
xlabel('Iterations');
legend('EXP3','best arm','UCB','OMA');
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Qu 2: EXP3 vs. tuned epsilon-Greedy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. By increasing T=4*10^4 or higher, we can see that EXP3 will beat 
% it eventually, T>>1, given that cumulative regret of EXP3 is sublinear, 
% whereas this is not the case for epsilon-Greedy
% 2. The slowness of EXP3 is explained by a slow regret decay adapted to the
% adversarial worst-case scenario but too slow for the stochastic setting
% 3. Nevertheless, in practice for finite horizons, epsilon-Greedy can beat
% EXP3, if epsilon can be finely tuned and in stochastic settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

