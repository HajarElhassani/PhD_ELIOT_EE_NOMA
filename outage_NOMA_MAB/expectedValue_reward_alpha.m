%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the expected reward for a given arm and a given alpha 'exhaustive search'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: expectedValue -> expected reward for a given arm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: selected_arm -> input arm; 
%         Pmax -> available transmit power at the base station
%         threshold 1 and 2 -> the SNR minimum thresholds
%         sigma1, sigma2, var_h1, var_h2 -> noise and channel links variances
%         of user 1 and user 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function expectedValue = expectedValue_reward_alpha(user,alpha,threshold1,threshold2,Pmax,sigma1,sigma2,var_h1,var_h2)

% User 1 performs SIC decoding and User 2 does SUD
if user==1 
    P=[Pmax*alpha Pmax*(1-alpha)];
    value1=sigma1*threshold1/P(1);
    value2=sigma1*threshold2/(P(2)-threshold2*P(1));
    value3=sigma2*threshold2/(P(2)-threshold2*P(1));
    variance1=var_h1;
    variance2=var_h2;
else
    % User 2 performs SIC decoding and User 1 does SUD
    P=[Pmax*(1-alpha) Pmax*alpha];
    value1=sigma2*threshold2/P(2);
    value2=sigma2*threshold1/(P(1)-threshold1*P(2));
    value3=sigma1*threshold1/(P(1)-threshold1*P(2));
    variance1=var_h2;
    variance2=var_h1;
end
% expected reward from the analytical expressions
if value2<0
    expectedValue=0;
else
    expectedValue=exp(-max(value1,value2)/(2*variance1))*exp(-value3/(2*variance2));
end

end 