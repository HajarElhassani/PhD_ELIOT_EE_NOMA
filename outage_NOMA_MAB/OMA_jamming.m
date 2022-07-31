function [outage,R_cum] = OMA_jamming(t,Pmax,G,threshold1,threshold2,R_cum,jam_user1,jam_user2)

%calcul gu reward
gain=reward_OMA_jamming(Pmax,G,threshold1,threshold2,jam_user1,jam_user2);
%calcul du reward cum
R_cum=R_cum+gain;
%calcul du outage
outage=1-R_cum/t;

end

