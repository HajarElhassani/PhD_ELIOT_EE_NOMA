%% generate channels between BS-users links
%input:      order -> (order=1: if SIC ordering should be performed;
%                     order=0: if SIC ordering should not be performed
%            cor_user -> coordinates of users
%            alpha -> pathlass exponent
%            sigma -> noise variance
%            error -> channel error (=0 in the case of perfect channel)
%output:     H= channel 
%            I -> the SIC ordering

function [H,I]=channelGain_BS(order,cor_user,alpha,sigma,error)

I=[];
%% Distances
%Distance between users and BSs
Dis_user_BS=abs(cor_user); 
h=(sqrt(Dis_user_BS.^(-alpha))-error)./sqrt(sigma);
if order==1
    [H,I]=sort(h,'descend');
else
    H=h;
end

