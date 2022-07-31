%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the best arm that maximizes UCB of the empirical_reward at time t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: arm_ind -> index of optimal arm
%         arm -> optimal arm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: t -> current iteration index;
%         alpha -> UCB parameter exploration vs. exploitation
%         arms -> set of arms; 
%         number_of_selected_arms -> ax1 recursive vector containing the 
%         number of times each arm has been selected at the beginning of
%         iteration t;
%         empirical_mean -> ax1 recursive vector containing the empirical mean 
%         rewards of each arm at the beginning of iteration t;
%         R_cum -> recursive scalar of updated cumulated rewards at the beginning 
%         of iteration t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [arm,arm_ind] = argmax_UCB(t, alpha, arms, empirical_mean, number_of_selected_arms)

    % Upper confidence bound at time t for all arms
    Q = empirical_mean + sqrt(alpha*log(t)./(2*number_of_selected_arms));

    % find the best arm maximizing Q
    [val,idx]=find(Q==max(Q));
    m=ceil(rand*length(idx));
    arm_ind=idx(m);
    arm=arms(arm_ind,:);

end

