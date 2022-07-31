%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration t of UCB algorithm, for t>a, after all arms have been tried
% once
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: cum_reward -> cummulative reward for the horizon 1 : t;
%         EE -> energy efficiency for the horizon 1 : t;
%         number_of_selected_arms -> ax1 recursive vector of the updated 
%         number of times each arm has been selected at the end of
%         iteration t;
%         empirical_mean -> ax1 recursive vector containing the updated 
%         empirical mean rewards of each arm at the end of iteration t;
%         R_cum -> recursive scalar of updated cumulated rewards at the end 
%         of iteration t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: t -> current iteration index;
%         arms -> set of arms; 
%         alpha -> UCB parameter exploration vs. exploitation
%         G = G(t,:) -> a vector of size 1x2 of the channels of thw two 
%         users at time t;
%         Pmax -> available transmit power at the base station
%         Pc -> circuit power
%         threshold 1 and 2 -> the SNR minimum thresholds
%         number_of_selected_arms -> ax1 recursive vector containing the 
%         number of times each arm has been selected at the beginning of
%         iteration t;
%         empirical_mean -> ax1 recursive vector containing the empirical mean 
%         rewards of each arm at the beginning of iteration t;
%         R_cum -> recursive scalar of updated cumulated rewards at the beginning 
%         of iteration t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cum_reward,EE,number_of_selected_arms,empirical_mean,R_cum] = UCB(t,arms,alpha,Pmax,Pc,G,threshold1,threshold2,number_of_selected_arms,empirical_mean,R_cum)

% optimal arm maximizing the upper confidence bound of the empirical average reward
delta=upper_bound(alpha,t,number_of_selected_arms);
[arm,arm_ind]=argmax(arms,empirical_mean,delta);

%get the feedback from the chosen arm
feedback=reward(arm,G,Pmax,threshold1,threshold2);
feedback_EE=(log2(1+threshold1)+log2(1+threshold2))*feedback/(Pmax*arms(arm_ind,2)+Pc);

%compute thge cumulative reward of the chosen arm
R_cum=R_cum+feedback_EE;
cum_reward=R_cum;

% update the number of times an arm is chosen and the empirical mean of the chosen arm
number_of_selected_arms(arm_ind)=number_of_selected_arms(arm_ind)+1;
empirical_mean(arm_ind)=empirical_mean_update(empirical_mean(arm_ind),feedback_EE,number_of_selected_arms(arm_ind));

%energy efficiency
EE=R_cum/t;

end