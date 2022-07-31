%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouput: feedback_chosen -> 1-bit feedback or instantaneous reward as a result playing the chosen arm
%        feedback_not_chosen -> 1-bit feedback or instantaneous reward as a result playing the non-chosen arm   
%       jam_user1,jam_user2 -> SNR of the jammer on respectively the user 1 and user2
% Inputs: arm played
%       G = G(t,:) -> a vector of size 1x2 with the channels of thw two users at time t
%       Pmax - available transmit power at the base station
%       threshold 1 and 2 -> the SNR minimum thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [feedback_chosen,feedback_not_chosen,jam_user1,jam_user2] = reward_jamming(arm,G,Pmax,threshold1,threshold2)

if arm(1)==1
    P=[Pmax*arm(2) Pmax*(1-arm(2))];
    %recalculate SNR (action 1)
    SNR1_arm1=G(1)*P(1);
    SNR2_arm1=G(2)*P(2)/(G(2)*P(1)+1);
    SNR_SIC_arm1=G(1)*P(2)/(G(1)*P(1)+1);
    %recalculate SNR (action 2)
    SNR1_arm2=G(1)*P(2)/(G(1)*P(1)+1);
    SNR2_arm2=G(2)*P(1);
    SNR_SIC_arm2=G(2)*P(2)/(G(2)*P(1)+1);
    
    %recalculate the reward (action 1)
    if (SNR1_arm1<=threshold1) || (min(SNR2_arm1,SNR_SIC_arm1)<=threshold2)
        feedback_chosen=0;
    else
        feedback_chosen=1;
    end
    
    %recalculate the reward (action 2)
    if (SNR2_arm2<=threshold2) || (min(SNR1_arm2,SNR_SIC_arm2)<=threshold1)
        feedback_not_chosen=0;
    else
        feedback_not_chosen=1;
    end
    
else
    P=[Pmax*(1-arm(2)) Pmax*arm(2)];
    %recalculate SNR (action 2)
    SNR1_arm2=G(1)*P(1)/(G(1)*P(2)+1);
    SNR2_arm2=G(2)*P(2);
    SNR_SIC_arm2=G(2)*P(1)/(G(2)*P(2)+1);
    
    %recalculate SNR (action 1)
    SNR1_arm1=G(1)*P(2);
    SNR2_arm1=G(2)*P(1)/(G(2)*P(2)+1);
    SNR_SIC_arm1=G(1)*P(1)/(G(1)*P(2)+1);
    
    %recalculate the reward (action 2)
    if (SNR2_arm2<=threshold2) || (min(SNR1_arm2,SNR_SIC_arm2)<=threshold1)
        feedback_chosen=0;
    else
        feedback_chosen=1;
    end
    
    %recalculate the reward (action 1)
    if (SNR1_arm1<=threshold1) || (min(SNR2_arm1,SNR_SIC_arm1)<=threshold2)
        feedback_not_chosen=0;
    else
        feedback_chosen=1;
    end
end

%% Since the reward of the chosen arm is 0, the jammer doesn't have to attack to kill the reward of arm
if (feedback_chosen==0)&&(feedback_not_chosen==0)
    jam_user1=0;
    jam_user2=0;
elseif (feedback_chosen==0)&&(feedback_not_chosen==1)
    jam_user1=0;
    jam_user2=0;
    
%% The jammer attacks to kill the reward of arm
elseif(feedback_chosen==1)&&(feedback_not_chosen==0)
    [feedback_chosen,feedback_not_chosen,jam_user1,jam_user2]=jamming_SNR(arm,Pmax,G,threshold1,threshold2);
elseif(feedback_chosen==1)&&(feedback_not_chosen==1)
    [feedback_chosen,feedback_not_chosen,jam_user1,jam_user2]=jamming_SNR(arm,Pmax,G,threshold1,threshold2);
end


end