%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouput: 1-bit feedback or instantaneous reward as a result of playing a certain arm 
% Inputs: arm played
%       Gt = G(t,:) -> a vector of size 1x2 with the channels of thw two users at time t
%       Pmax - available transmit power at the base station
%       threshold 1 and 2 -> the SNR minimum thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [feedback] = reward(arm,Gt,Pmax,threshold1,threshold2)


if arm(1)==1     

    % User 1 does SIC decoding, User 2 does SUD
    P=[Pmax*arm(2) Pmax*(1-arm(2))];
    
    SNR1=Gt(1)*P(1);
    SNR2=Gt(2)*P(2)/(Gt(2)*P(1)+1);
    SNR_SIC=Gt(1)*P(2)/(Gt(1)*P(1)+1);
    
    % Compute the bit of feedback
    if (SNR1<threshold1) || (min(SNR2,SNR_SIC)<threshold2)
        feedback=0; % if in outage feedback = 0
    else
        feedback=1;
    end

else
    
    % User 2 does the SIC decoding, User 1 does SUD 
    P=[Pmax*(1-arm(2)) Pmax*arm(2)];
    
    SNR1=Gt(1)*P(1)/(Gt(1)*P(2)+1);
    SNR2=Gt(2)*P(2);
    SNR_SIC=Gt(2)*P(1)/(Gt(2)*P(2)+1);
    
    % Compute the bit of feedback
    if (SNR2<threshold2) || (min(SNR1,SNR_SIC)<threshold1)
        feedback=0; % if in outage feedback = 0
    else
        feedback=1;
    end
end

end