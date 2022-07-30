%% This function computes the optimal reflection coefficient in the case of four backscatter device
%input: channels
%output: rho:optimal reflection coefficients vector and the corresponding backscattered+direct channel
function [rho,H_NOMA] = optimal_rho_4BD(G_users,G_BS_BD1,G_BS_BD2,G_BS_BD3,G_BS_BD4,G_BD1_users,G_BD2_users,G_BD3_users,G_BD4_users)

%% find the optimal rho1, rho2, rho3, rho4 using EXHAUSTIVE SEARCH  

%channel of the backscattered link from backscatter 1,2,3,4
G_BD1=G_BS_BD1*G_BD1_users;
G_BD2=G_BS_BD2*G_BD2_users;
G_BD3=G_BS_BD3*G_BD3_users;
G_BD4=G_BS_BD4*G_BD4_users;

[Rho1,Rho2,Rho3,Rho4]=ndgrid(0:.01:1);
rho1=Rho1(:);
rho2=Rho2(:);
rho3=Rho3(:);
rho4=Rho4(:);

%take the maximum value of rho1, rho2, rho3 that satisfy the constraints
solidx = find(all([0 <= rho1,0 <= rho2, 0 <= rho3, 0 <= rho4, rho1 <= 1,rho2 <= 1, rho3 <= 1,rho4 <= 1,((G_users(1)-G_users(2))+(G_BD1(1)-G_BD1(2))*sqrt(rho1)+(G_BD2(1)-G_BD2(2))*sqrt(rho2)+(G_BD3(1)-G_BD3(2))*sqrt(rho3)+(G_BD4(1)-G_BD4(2))*sqrt(rho4))>=0], 2));
rho=[max(rho1(solidx)),max(rho2(solidx)),max(rho3(solidx)),max(rho4(solidx))];

%channel gains for the backscattered link --> (backscatter 1,2,3,4) under OMA
H_NOMA=(G_users+sqrt(rho(1))*G_BD1+sqrt(rho(2))*G_BD2+sqrt(rho(3))*G_BD3+sqrt(rho(4))*G_BD4).^2;

end

