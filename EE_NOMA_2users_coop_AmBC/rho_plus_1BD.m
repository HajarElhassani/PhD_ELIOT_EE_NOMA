%% This function computes the optimal reflection coefficient for one backscatter device
%% input: channel gains
%% output: the reflection coefficient rho

function rho_vect = rho_plus_1BD(G_users,G_BS_BD,G_BD_users)

% the backscattered channel
G_bar=G_BS_BD.*G_BD_users;

% return empty if the condition G_bar(k)-G_bar(k-1)>0 is not satisfied 
rho_vect=[];

for k=2:length(G_BD_users)
    if (G_bar(k)>G_bar(k-1))
        rho_vect(k-1)=(G_users(k-1)-G_users(k))/(G_bar(k)-G_bar(k-1));
    end
end

end

