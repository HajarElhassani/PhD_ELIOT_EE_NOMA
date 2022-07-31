%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: cum_reward -> cummulative reward for the horizon 1 : t;
%         EE -> outage for the horizon 1 : t;
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
%         G = G(t,:) -> a vector of size 1x2 of the channels of thw two 
%         users at time t;
%         Pmax -> available transmit power at the base station
%         Pc -> circuit power
%         threshold 1 and 2 -> the SNR minimum thresholds
%         weights -> ax1 recursive mutliplicative weights of EXP3 at the
%         begining of iteration t
%         number_of_selected_arms -> ax1 recursive vector containing the 
%         number of times each arm has been selected at the beginning of
%         iteration t;
%         R_cum -> recursive scalar of updated cumulated rewards at the beginning 
%         of iteration t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [cum_reward,EE,weights,number_of_selected_arms,empirical_mean,R_cum] = EXP3(t,arms,Pmax,Pc,weights,gamma,eta,G,threshold1,threshold2,number_of_selected_arms,empirical_mean,R_cum)

% update the probability distribution of the arms folowing exp3 
probability_distribution=distr(weights,gamma);

% draw a random arm following 
arm_ind=gendist(probability_distribution,1,1);
arm=arms(arm_ind,:);

%feedback of the chosen arm
feedback=reward(arm,G,Pmax,threshold1,threshold2);
feedback_EE=(log2(1+threshold1)+log2(1+threshold2))*feedback/(Pmax*arms(arm_ind,2)+Pc);

%cumulated reward
R_cum=R_cum+feedback_EE;
cum_reward=R_cum;

% unbiased estimation of the estimated arm
estimate_gain=feedback_EE/probability_distribution(arm_ind);

% multiplicative weights update
weights(arm_ind)=weights(arm_ind)*exp(estimate_gain*eta);

%energy efficiency
EE=R_cum/t;

end