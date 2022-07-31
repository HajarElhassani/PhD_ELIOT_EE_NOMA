%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update of empirical mean 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouput: empirical_mean_new -> scalar, update of empirical average reward 
% of a chosen arm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%         empirical_mean -> scalar, previous value of empirical average
%         reward of the arm
%         feedback -> one bit reward as a result of choosing the arm
%         Na -> scalar, number of times the arm has been chosen in total
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function empirical_mean = empirical_mean_update(mean,gain,Na)

    empirical_mean=(1-1/Na)*mean+gain/Na;

end

