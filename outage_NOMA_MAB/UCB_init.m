
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration t of UCB algorithm, for t>a, after all arms have been tried
% once
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: cum_reward_t -> cummulative reward for the horizon 1 : t;
%         outage_t -> outage for the horizon 1 : t;
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
%         Gt = G(t,:) -> a vector of size 1x2 of the channels of thw two 
%         users at time t;
%         Pmax -> available transmit power at the base station
%         threshold 1 and 2 -> the SNR minimum thresholds
%         vepsilon -> parameter of epsilon-Greedy, probability of
%         exploration;
%         number_of_selected_arms -> ax1 recursive vector containing the 
%         number of times each arm has been selected at the beginning of
%         iteration t;
%         empirical_mean -> ax1 recursive vector containing the empirical mean 
%         rewards of each arm at the beginning of iteration t;
%         R_cum -> recursive scalar of updated cumulated rewards at the beginning 
%         of iteration t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [cum_reward_t,outage_t,number_of_selected_arms,empirical_mean,R_cum] = UCB_init(t,arms,Pmax,Gt,threshold1,threshold2,number_of_selected_arms,empirical_mean,R_cum)

    % choose the arm
    arm_ind = t;   % we start by playing arms t for t = 1 to a 
    arm=arms(arm_ind,:); 

    %get the feedback from the chosen arm
    feedback = reward(arm,Gt,Pmax,threshold1,threshold2);

    %compute thge cumulative reward of the chosen arm
    R_cum=R_cum+feedback;
    cum_reward_t=R_cum;

    % update the number of times an arm is chosen and the empirical mean of the chosen arm
    number_of_selected_arms(arm_ind)=number_of_selected_arms(arm_ind)+1;
    empirical_mean(arm_ind)=empirical_mean_update(empirical_mean(arm_ind),feedback,number_of_selected_arms(arm_ind));

    %outage
    outage_t=1-R_cum/t;

end

