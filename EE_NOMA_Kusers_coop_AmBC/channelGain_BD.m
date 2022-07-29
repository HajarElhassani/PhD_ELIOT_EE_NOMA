%% This function generates the channels between the bacackscatter and users

function g=channelGain_BD(cor_BD,cor_user,alpha)
Dis_BD_users=abs((real(cor_user)-real(cor_BD))+i*(imag(cor_user)-imag(cor_BD))); 
g=sqrt(Dis_BD_users.^(-alpha));
end
