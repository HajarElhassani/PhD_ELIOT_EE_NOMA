%% This function computes the optimal reflection coefficient in the case of one backscatter device
%input: channels
%output: rho:optimal reflection coefficient and the corresponding backscattered+direct channel

function [rho_NOMA,H_NOMA] = optimal_rho_1BD(G_users,G_BS_BD,G_BD_users)
%the channel of the backscattered link
G_bar=G_BS_BD*G_BD_users;

%returns empty if the condition G_bar(k)-G_bar(k-1)>0 is not satisfied
rho=[];

%% compute the optimal rho 
if (G_bar(2)>G_bar(1))
    rho=((G_users(1)-G_users(2))/(G_bar(2)-G_bar(1)))^2;
end


if (isempty(rho))
    rho_NOMA=1;
else
    rho_NOMA=min(1,rho);
end

%channel gain of the backscattered link (one backscatter) --> reflection coefficient computed for NOMA
H_NOMA=(G_users+sqrt(rho_NOMA)*G_bar).^2;

end

