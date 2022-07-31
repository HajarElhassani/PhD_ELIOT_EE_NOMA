%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the best expected reward, and of the worst expected reward
% for certain alpha(fraction of the power) 'exhaustive search'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: mu_best -> best expected reward over the arms
%         mu_worst -> worst expected reward over the arms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: users -> the set of users (who performs SIC); 
%         alpha -> the fraction of power 
%         Pmax -> available transmit power at the base station
%         threshold 1 and 2 -> the SNR minimum thresholds
%         sigma1, sigma2, var_h1, var_h2 -> noise and channel links variances
%         of user 1 and user 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [mu_best,mu_worst] = expected_offline_policy_alpha(users,alpha,Pmax,threshold1,threshold2,sigma1,sigma2,var_h1,var_h2)

% Expected reward calculation for each arm
for i=1:length(users)
mu_a(i)=expectedValue_reward_alpha(users(i),alpha,threshold1,threshold2,Pmax,sigma1,sigma2,var_h1,var_h2);
end

%choose the arm with the maximum statistical mean of reward
[val1,idx1]=find(mu_a==max(mu_a));
m1=ceil(rand*length(idx1)); %Choose a max at random in case we have more than 2 arms
ind1=idx1(m1);

%choose the arm with the minimum statistical mean of reward
[val2,idx2]=find(mu_a==min(mu_a));
m2=ceil(rand*length(idx2)); %Choose a min at random in case we have more than 2 arms
ind2=idx2(m2);

mu_best=mu_a(ind1);
mu_worst=mu_a(ind2);

end

