%% This function computes the optimal reflection coefficient in the case of two backscatter device
%input: channels
%output: rho:optimal reflection coefficients vector and the corresponding backscattered+direct channel

function [rho,H_NOMA,H_OMA] = optimal_rho_2BD(G_users,G_BS_BD1,G_BS_BD2,G_BD1_users,G_BD2_users)

%channel of the backscattered link from backscatter 1 and 2
G1=G_BS_BD1*G_BD1_users;
G2=G_BS_BD2*G_BD2_users;

%% different solutions according to the 4 cases in the theorem 
rho1_bar=(G_users(1)-G_users(2))/(G1(2)-G1(1))^2;
rho2_bar=(G_users(1)-G_users(2))/(G2(2)-G2(1))^2;

rho1_tild=(((G_users(1)-G_users(2))-(G2(2)-G2(1)))/(G1(2)-G1(1)))^2;
rho2_tild=(((G_users(1)-G_users(2))-(G1(2)-G1(1)))/(G2(2)-G2(1)))^2;


%% discuss according to the channel gain cases
if ((G1(2)-G1(1))<=0 && (G2(2)-G2(1))<=0)
    rho_1=1;
    rho_2=1;
elseif ((G1(2)-G1(1))<=0 && (G2(2)-G2(1))>0)
    rho_1=1;
    rho_2=min(1,rho2_tild);
elseif ((G1(2)-G1(1))>0 && (G2(2)-G2(1))<=0)
    rho_1=min(1,rho1_tild);
    rho_2=1;
else
    if (G2(2)*G1(1)-G1(2)*G2(1))>=0
        rho_1=min(1,rho1_bar);
        rho_2=min(1,max(0,rho2_tild));
    else
        rho_1=min(1,max(0,rho1_tild));
        rho_2=min(1,rho2_bar);
    end
end

%optimal reflection coefficient vector (backscatter 1, backscatter 2)
rho=[rho_1,rho_2];

%channel gains for the backscattered link --> (backscatter 1, backscatter2) under NOMA
H_NOMA=(G_users+sqrt(rho(1))*G1+sqrt(rho(2))*G2).^2;

%channel gains for the backscattered link --> (backscatter 1, backscatter2) under OMA
H_OMA=(G_users+G1+G2).^2;

end