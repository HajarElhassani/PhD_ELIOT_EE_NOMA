%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the expected reward for OMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: mu -> expected reward for OMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%         Pmax -> available transmit power at the base station
%         threshold 1 and 2 -> the SNR minimum thresholds
%         sigma1, sigma2, var_h1, var_h2 -> noise and channel links variances
%         of user 1 and user 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mu = expected_value_OMA(threshold1,threshold2,Pmax,sigma1,sigma2,var_h1,var_h2)

    mu = exp(-sigma1*threshold1/(2*Pmax*var_h1))*exp(-sigma2*threshold2/(2*Pmax*var_h2));

end

