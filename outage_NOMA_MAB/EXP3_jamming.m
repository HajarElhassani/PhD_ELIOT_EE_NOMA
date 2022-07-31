%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: 
%         outage -> outage for the horizon 1 : t;
%         weights -> ax1 recursive vector of the updated multiplicative 
%         weights of EXP3 at the end of iteration t
%         R_cum -> recursive scalar of updated cumulated rewards at the end 
%         of iteration t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: t -> current iteration index;
%         gamma,eta -> EXP3 parameters exploration vs. exploitation
%         weights -> ax1 recursive mutliplicative weights of EXP3 at the
%         begining of iteration t
%         R_cum -> recursive scalar of updated cumulated rewards at the beginning 
%         of iteration t
%         feedback_vect -> feedback vector for the chosen and non-chosen arm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cum_feedback,outage,weights,R_cum] = EXP3_jamming(t,eta,weights,gamma,feedback_vect,R_cum)

% update the probability distribution of the arms folowing exp3 
probability_distribution=(1-gamma)*(weights/sum(weights))+gamma/length(weights);
% draw a random arm following 
arm_ind=((find(rand<cumsum(probability_distribution),1,'first')));

%feedback of the chosen arm
if arm_ind==1
    feedback=feedback_vect(1);
else
    feedback=feedback_vect(2);
end

%cumulated reward
R_cum=R_cum+feedback;
cum_feedback=R_cum;

% unbiased estimation of the estimated arm
estimate_gain=feedback/probability_distribution(arm_ind);

% multiplicative weights update
weights(arm_ind)=weights(arm_ind)*exp(estimate_gain*eta);
%weights(arm_ind)=weights(arm_ind)*exp(estimate_gain*gamma/length(weights));

%outage
outage=1-R_cum/t;

end