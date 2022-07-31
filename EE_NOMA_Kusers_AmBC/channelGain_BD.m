%% generate channels between BD-users links
%input:      cor_BD -> coordinates of the backscatter device
%            cor_user -> coordinates of the users
%            alpha -> pathlass exponent
%            error -> channel error (=0 in the case of perfect channel)
%output:     g= channel 

function g=channelGain_BD(cor_BD,cor_user,alpha,error)
%% Distances
%Distance between users and BSs
Dis_BD_users=abs((real(cor_user)-real(cor_BD))+i*(imag(cor_user)-imag(cor_BD)));
g=sqrt(Dis_BD_users.^(-alpha))-error;
end
