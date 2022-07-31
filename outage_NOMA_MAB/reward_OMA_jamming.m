function feedback = reward_OMA_jamming(Pmax,G,threshold1,threshold2,jam_user1,jam_user2)

SNR1=G(1)*Pmax/(1+jam_user1);
SNR2=G(2)*Pmax/(1+jam_user2);

if (SNR1<threshold1) || (SNR2<threshold2)
    feedback=0;
else
    feedback=1;
end
    
end

