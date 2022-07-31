%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouput: feedback_chosen -> 1-bit feedback or instantaneous reward as a result playing the chosen arm
%        feedback_not_chosen -> 1-bit feedback or instantaneous reward as a result playing the non-chosen arm   
%       jam_user1,jam_user2 -> SNR of the jammer (noise introduced by the jammer) on respectively the user 1 and user2
% Inputs: arm played
%       G = G(t,:) -> a vector of size 1x2 with the channels of thw two users at time t
%       Pmax - available transmit power at the base station
%       threshold 1 and 2 -> the SNR minimum thresholds
%       arm by the jammer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [feedback_chosen,feedback_not_chosen,jam_user1,jam_user2] = jamming_SNR(arm,Pmax,G,threshold1,threshold2)

%action (a1) selected: user 1 performs SIC
if(arm(1)==1)
    P=[Pmax*arm(2) Pmax*(1-arm(2))];
    % compute the SNR of the jammer on user 1 and user 2
    jam_user1=min(P(2)*G(1)/threshold2-(1+P(1)*G(1)),P(1)*G(1)/threshold1-1);
    jam_user2=P(2)*G(2)/threshold2-(1+P(1)*G(2));
    if (jam_user1 <=0 || jam_user2<= 0)
        error  ('No negative numbers allowed');
    end
    
    %recalculate SNR (action 1)
    SNR1_arm1=G(1)*P(1)/(1+jam_user1);
    SNR2_arm1=G(2)*P(2)/(G(2)*P(1)+jam_user2+1);
    SNR_SIC_arm1=G(1)*P(2)/(G(1)*P(1)+jam_user1+1);
    %recalculate the reward (action 1)
    if (SNR1_arm1<=threshold1) || (min(SNR2_arm1,SNR_SIC_arm1)<=threshold2)
        feedback_chosen=0;
    else
        feedback_chosen=1;
    end
    
    %recalculate SNR (action 2)
    SNR1_arm2=G(1)*P(2)/(G(1)*P(1)+jam_user1+1);
    SNR2_arm2=G(2)*P(1)/(1+jam_user2);
    SNR_SIC_arm2=G(2)*P(2)/(G(2)*P(1)+jam_user2+1);
    %recalculate the reward (action 2)
    if (SNR2_arm2<=threshold2) || (min(SNR1_arm2,SNR_SIC_arm2)<=threshold1)
        feedback_not_chosen=0;
    else
        feedback_not_chosen=1;
    end
    
    
    
      
    
% action (2) selected: user 2 performs SIC    
else
    P=[Pmax*(1-arm(2)) Pmax*arm(2)];
    % compute the SNR of the jammer on user 1 and user 2 SNR of jammer_user1
    jam_user1=P(1)*G(1)/threshold1-(1+P(2)*G(1));
    jam_user2=min(P(2)*G(2)/threshold2-1,P(1)*G(2)/threshold1-(1+P(2)*G(2)));
    
    if (jam_user1 <=0 || jam_user2<= 0)
        error  ('No negative numbers allowed');
    end
    
    
    %recalculate SNR (action 2)
    SNR1_arm2=G(1)*P(1)/(G(1)*P(2)+jam_user1+1);
    SNR2_arm2=G(2)*P(2)/(1+jam_user2);
    SNR_SIC_arm2=G(2)*P(1)/(G(2)*P(2)+jam_user2+1);
    
    %recalculate SNR (action 1)
    SNR1_arm1=G(1)*P(2)/(1+jam_user1);
    SNR2_arm1=G(2)*P(1)/(G(2)*P(2)+jam_user2+1);
    SNR_SIC_arm1=G(1)*P(1)/(G(1)*P(2)+jam_user1+1);
    
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
        feedback_not_chosen=1;
    end
   
    
end


end


