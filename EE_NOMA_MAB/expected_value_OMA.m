
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the expected reward for OMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: mu -> expected reward for OMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: 
%         selected_arm -> input arm
%         Pmax -> available transmit power at the base station
%         threshold 1 and 2 -> the SNR minimum thresholds
%         sigma1, sigma2, var_h1, var_h2 -> noise and channel links variances
%         of user 1 and user 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mu = expected_value_OMA(selected_arm,threshold1,threshold2,Pmax,sigma1,sigma2,variance1,variance2)

P=[Pmax*selected_arm(2) Pmax*selected_arm(2)];
value1=sigma1*threshold1/P(1);
value2=sigma2*threshold2/P(2);


mu=exp(-value1/(2*variance1))*exp(-value2/(2*variance2));

end





