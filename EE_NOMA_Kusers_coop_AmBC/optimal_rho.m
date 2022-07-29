%% This function compute Pmin in NOMA with BD and OMA with BD + computing optimal rho

%% input: channel gains, A=2^(2*Rmin), R
%% output: optimal rho (rho_NOMA), Pmin for NOMA (Pmin_NOMA),Pmin for OMA (Pmin_OMA), Gamma_k according to notations in the paper (G_NOMA, G_OMA)

function [rho_NOMA,G_NOMA,Pmin_NOMA,G_OMA,Pmin_OMA] = optimal_rho(G_BS_users,G_BS_BD,G_BD_users,A,R)
%% get the optimal tho according to theorem 1
if (isempty(R))
    rho_NOMA=1;
else
    rho_NOMA=min(1,min(R));
end

%% compute the channel gain of the backscattered link 
G_NOMA=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2;
G_OMA=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2;

%% compute Pmin
Pmin_NOMA=0;
for i=1:length(G_BS_users)
    Pmin_NOMA=Pmin_NOMA+(A(i)-1)/G_NOMA(i)*prod(A(i+1:length(G_BS_users)));
end
Pmin_OMA=sum((A.^length(G_BS_users)-1)./G_OMA);

end
