
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the expected reward for OMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: mu_best-> expected reward for OMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%         arms -> input arms
%         Pmax -> available transmit power at the base station
%         threshold 1 and 2 -> the SNR minimum thresholds
%         sigma1, sigma2, var_h1, var_h2 -> noise and channel links variances
%         of user 1 and user 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mu_best = expected_offline_policy_OMA(arms,Pmax,Pc,threshold1,threshold2,sigma1,sigma2,var_h1,var_h2)

%compute the reward of each arm
for i=1:length(arms)
mu(i)=expected_value_OMA(arms(i,:),threshold1,threshold2,Pmax,sigma1,sigma2,var_h1,var_h2);
end
mu_a=(1/2*log2(1+threshold1)+1/2*log2(1+threshold2))*mu./(Pmax*arms(:,2)'+Pc);

%choose the arm with the best statistical mean of reward
[val1,idx1]=find(mu_a==max(mu_a));
m1=ceil(rand*length(idx1)); %Choose a max at random in case we have more than 2 arms
ind1=idx1(m1);

mu_best=mu_a(ind1);
end

