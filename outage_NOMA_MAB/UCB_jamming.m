%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: same as UCB 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: same as UCB +
%         reward_vec -> reward of chosen and not-chosen arm; 
%         feedback_not_chosen -> feedback of the not-chosen arm
%         jam_user1,jam_user2 -> SNR of the jammer on respectively the user 1 and user2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cum_reward,outage,feedback_vec,feedback_not_chosen,number_of_selected_arms,empirical_mean,R_cum,jam_user1,jam_user2] = UCB_jamming(t,arms,alpha,Pmax,G,threshold1,threshold2,number_of_selected_arms,empirical_mean,R_cum)

% optimal arm maximizing the upper confidence bound of the empirical
% average reward
[arm, arm_ind] = argmax_UCB(t, alpha, arms,empirical_mean,number_of_selected_arms);

% compute the rewards of the chosen and non-chosen arm and SNR of jammer
[feedback_chosen,feedback_not_chosen,jam_user1,jam_user2]=reward_jamming(arm,G,Pmax,threshold1,threshold2);

%save the rewards in order to have (feedback(a1 -> chosen arm),feedback(a2 -> non-chosen arm))
if arm_ind==1
    feedback_vec=[feedback_chosen,feedback_not_chosen];
else
    feedback_vec=[feedback_not_chosen,feedback_chosen];
end

%compute the cumulative reward of the chosen arm
R_cum=R_cum+feedback_chosen;
cum_reward=R_cum;

% update the number of times an arm is chosen and the empirical mean of the chosen arm
number_of_selected_arms(arm_ind)=number_of_selected_arms(arm_ind)+1;
empirical_mean(arm_ind)=empirical_mean_update(empirical_mean(arm_ind),feedback_chosen,number_of_selected_arms(arm_ind));

%compute  the outage
outage=1-R_cum/t;


end

