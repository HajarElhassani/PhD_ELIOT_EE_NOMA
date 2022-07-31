clear all;
close all;
clc;

%% Parameter initialization
N=1e3; % number of parallel runs for the average regret
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

% arms
arms=[1 0.4;2 0.4];
a=length(arms); % number of arms

% EXP3 learning parameters
gamma=min(1,sqrt(a*log(a)/((exp(1)-1)*T)));
eta=0.02; % gamma/a;

% Parameters of learning algorithms 
alpha=1; % upper bound parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iteration for N runs 
for n=1:N
    % generating channel gains for each run 
    [G] = channel_gain(T,sigma1,sigma2,var_h1,var_h2);  
    %----------------------------------------------------------------------
    %% Initialization of learning algorithms
    %----------------------------------------------------------------------
    
    % Initialization for UCB algorithm
    number_of_sel_UCB=ones(1,a); % initial number of times arms were selected
    empirical_mean_UCB=zeros(1,a);% initial rewards sample (empirical reward) from each arm
    R_cum_UCB=0;
   
    
    % Initialisation for EXP3 algorithm
    number_of_sel_EXP3=ones(1,a);
    weights=ones(1,a);
    R_cum_EXP3=0;
    
    %OMA initialisation
    R_cum_OMA=0;    
    
    %% epsilon-greedy, UCB
    % Initialization
    % test each arm of index t for t = 1 : a
    %----------------------------------------------------------------------
    for t = 1 : a

        % call one iteration of UCB init
        [cum_reward_UCB(t),outage_UCB(t),reward_vect(t,:),reward_not_chosen(t),number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB,jam_user1,jam_user2]=UCB_init_jamming(t,arms,Pmax,G(t,:),threshold1_NOMA,threshold2_NOMA,number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB);
        [cum_reward_EXP3(t),outage_EXP3(t),weights,R_cum_EXP3] = EXP3_jamming(t,eta,weights,gamma,reward_vect(t,:),R_cum_EXP3);
        [outage_OMA(t),R_cum_OMA] = OMA_jamming(t,Pmax,G(t,:),threshold1_OMA,threshold2_OMA,R_cum_OMA,jam_user1,jam_user2);
        cum_reward_offline(t)=sum(reward_not_chosen(1:t));

    end
    %----------------------------------------------------------------------
    %% UCB
    % generic iterations for t > a
    %----------------------------------------------------------------------
    for t = a+1 : T
                  
        % call one generic iteration of UCB_
        [cum_reward_UCB(t),outage_UCB(t),reward_vect(t,:),reward_not_chosen(t),number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB,jam_user1,jam_user2]=UCB_jamming(t,arms,alpha,Pmax,G(t,:),threshold1_NOMA,threshold2_NOMA,number_of_sel_UCB,empirical_mean_UCB,R_cum_UCB);
        [cum_reward_EXP3(t),outage_EXP3(t),weights,R_cum_EXP3] = EXP3_jamming(t,eta,weights,gamma,reward_vect(t,:),R_cum_EXP3);
        [outage_OMA(t),R_cum_OMA] = OMA_jamming(t,Pmax,G(t,:),threshold1_OMA,threshold2_OMA,R_cum_OMA,jam_user1,jam_user2);
        cum_reward_offline(t)=sum(reward_not_chosen(1:t));

    end
    
    
    %cumulative regret
    cum_regret_UCB(n,:)=cum_reward_offline-cum_reward_UCB;
    cum_regret_EXP3(n,:)=cum_reward_offline-cum_reward_EXP3;
    
    %regret
    regret_UCB(n,:)=cum_regret_UCB(n,:)./(1:T);
    regret_EXP3(n,:)=cum_regret_EXP3(n,:)./(1:T);
    
    %outage 
    outage_UCB_n(n,:)=outage_UCB;
    outage_EXP3_n(n,:)=outage_EXP3;
    outage_off_n(n)=1-max(sum(reward_vect(:,1))/T,sum(reward_vect(:,2))/T);
    outage_OMA_n(n,:)=outage_OMA;


end

%% averaging over all runs 
mean_cum_regret_UCB=mean(cum_regret_UCB);
mean_cum_regret_EXP3=mean(cum_regret_EXP3);

mean_regret_UCB=mean(regret_UCB);
mean_regret_EXP3=mean(regret_EXP3);

mean_outage_UCB=mean(outage_UCB_n);
mean_outage_EXP3=mean(outage_EXP3_n);
mean_outage_offline=mean(outage_off_n);
mean_outage_OMA=mean(outage_OMA_n);


%% plotting figures

%cumulative regret
figure(1);
plot(1:T,mean_cum_regret_UCB,'o-', 'MarkerIndices', 1:1000:length(1:T));
hold on;
plot(1:T,mean_cum_regret_EXP3,'d-', 'MarkerIndices', 1:1000:length(1:T));
ylabel('Cumulative regret');
xlabel('Iterations');
legend('UCB','EXP3','location','best');
grid on;

%regret
figure(2);
plot(1:T,mean_regret_UCB,'o-', 'MarkerIndices', 1:1000:length(1:T));
hold on;
plot(1:T,mean_regret_EXP3,'d-', 'MarkerIndices', 1:1000:length(1:T));
ylabel('Regret');
xlabel('Iterations');
legend('UCB','EXP3','location','best');
grid on;

%outage
figure(3);
plot(1:T,100*mean_outage_EXP3,'-<','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 1:1000:length(1:T));
hold on;
plot(1:T,100*ones(1,T)*mean_outage_offline,'-*','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 1:1000:length(1:T));
hold on;
plot(1:T,100*mean_outage_UCB,'->', 'MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:1000:length(1:T));
hold on;
plot(1:T,100*mean_outage_OMA,'--k','MarkerSize',7,'LineWidth',1.3, 'MarkerIndices', 1:1000:length(1:T));
ylabel('Outage (%)');
xlabel('Iterations');
legend('EXP3','Best offline arm','UCB','OMA','location','best');
grid on;