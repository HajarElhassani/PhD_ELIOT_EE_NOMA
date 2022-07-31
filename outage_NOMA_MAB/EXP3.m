%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: cum_reward_t -> cummulative reward for the horizon 1 : t;
%         outage_t -> outage for the horizon 1 : t;
%         weights -> ax1 recursive vector of the updated multiplicative 
%         weights of EXP3 at the end of iteration t
%         number_of_selected_arms -> ax1 recursive vector of the updated 
%         number of times each arm has been selected at the end of
%         iteration t;
%         R_cum -> recursive scalar of updated cumulated rewards at the end 
%         of iteration t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: t -> current iteration index;
%         arms -> set of arms; 
%         gamma, eta -> EXP3 parameters exploration vs. exploitation
%         Gt = G(t,:) -> a vector of size 1x2 of the channels of thw two 
%         users at time t;
%         Pmax -> available transmit power at the base station
%         threshold 1 and 2 -> the SNR minimum thresholds
%         weights -> ax1 recursive mutliplicative weights of EXP3 at the
%         begining of iteration t
%         number_of_selected_arms -> ax1 recursive vector containing the 
%         number of times each arm has been selected at the beginning of
%         iteration t;
%         R_cum -> recursive scalar of updated cumulated rewards at the beginning 
%         of iteration t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [cum_reward_t, outage_t, weights,number_of_selected_arms,R_cum] = EXP3(t,gamma,eta,arms,Pmax,weights,Gt,threshold1,threshold2,number_of_selected_arms,R_cum)


% update the probability distribution of the arms folowing exp3 
probability_distribution=distr_exp3(weights,gamma);

% draw a random arm following 
arm_ind=gendist(probability_distribution,1,1);
arm=arms(arm_ind,:);

%feedback of the chosen arm
feedback=reward(arm,Gt,Pmax,threshold1,threshold2);

%cumulated reward
R_cum=R_cum+feedback;
cum_reward_t=R_cum;

% unbiased estimation of the estimated arm
estimate_gain=feedback/probability_distribution(arm_ind);

% multiplicative weights update
weights(arm_ind)=weights(arm_ind)*exp(estimate_gain*eta);

% outage
outage_t=1-R_cum/t;


end