clear all;
close all;
clc;

%% Input parameters

N=1e3;%number of channel realizations
K=3; %number of users
Pmax_dbm=80;%power budget at the BS in dbm
Pc_dbm=30;%circuit power in dbm
sigma_dbm=-20;%noise power in dbm

%call the function 'dbm_to_Watt' to convert from dbm to Watt
Pmax=dbm_to_Watt(Pmax_dbm);
Pc=dbm_to_Watt(Pc_dbm);
sigma=dbm_to_Watt(sigma_dbm);

%coordinates to create the cells
min_dis_BS_users=0.5;%the minimum distance between BS and users
radius_BS_users=20;%maximum distance between BS and users
min_dis_BS_BD=0.5;%the minimum distance between BD and BS
radius_BS_BD=4;%maximum distance between BD and BS
alpha=2.5;%pathloss exponent

%varying the minimum rate QoS constraint
Rmin=0:4;
A=(2.^(2*Rmin)).*ones(K,1);


%% compute the GEE as a function of Rmin


for n=1:N
    
    %% generating channels
    %generate x and y coordinates for users
    users_cordinates = coordinates(K,radius_BS_users,min_dis_BS_users)';
    %generate x and y coordinates for BD
    BD_cordinates = coordinates(1,radius_BS_BD,min_dis_BS_BD);
    
    %generate BS-BD channel
    G_BS_BD=channelGain_BS(BD_cordinates,alpha,sigma);
    %generate BS-users channels (in descending order -> SIC)
    [G_BS_users,I]=channelGain_BS(users_cordinates,alpha,sigma);
    G_BS_users=G_BS_users.^2;
    %generate BD-users channels before SIC order
    G_BD_users_unordered=channelGain_BD(BD_cordinates,users_cordinates,alpha);
    %ordering channels BD-users
    G_BD_users=G_BD_users_unordered(I);
    
    %% calling the function that computes R according to equation (6) in the paper
    R = rho_plus(sqrt(G_BS_users),G_BS_BD,G_BD_users);
    
    %% compute the minimum power budget (Pmin) required for meeting each QoS constraints in conventional NOMA and OMA
    for i=1:length(A)
        Pmin_OMA_conv(i)=sum((A(:,i).^K-1)./G_BS_users);
        
        Pmin_NOMA_conv(i)=(A(end,i)-1)/G_BS_users(end);
        for j=1:K
            Pmin_NOMA_conv(i)=Pmin_NOMA_conv(i)+(A(j,i)-1)/G_BS_users(j)*prod(A((j+1):K,i));
        end
        
        %% compute the optimal reflection coefficient(rho)
        if (isempty(R))
            rho_NOMA=1;
        else
            rho_NOMA=min(1,min(R));
        end
        
        %% compute Gamma according to the notations in the paper
        G_OMA_BD=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2;
        G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2;
        
        %% compute Pmin in OMA with BD and NOMA with BD
        Pmin_OMA_BD(i)=sum((A(i).^length(G_BS_users)-1)./G_OMA_BD);
        Pmin_NOMA_BD(i)=0;
        for j=1:length(G_BS_users)
            Pmin_NOMA_BD(i)=Pmin_NOMA_BD(i)+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:length(G_BS_users)));
        end
        
    end
    
    %% checking the feasability condition. If not satisfied regenerate channels
    while (min(Pmin_NOMA_conv)>Pmax || min(Pmin_OMA_conv)>Pmax || min(Pmin_NOMA_BD)>Pmax ||min(Pmin_OMA_BD)>Pmax)
        
        %generate x and y coordinates for users
        users_cordinates = coordinates(K,radius_BS_users,min_dis_BS_users)';
        %generate x and y coordinates for BD
        BD_cordinates = coordinates(1,radius_BS_BD,min_dis_BS_BD);
        
        %generate BS-BD channel
        G_BS_BD=channelGain_BS(BD_cordinates,alpha,sigma);
        %generate BS-users channels (in descending order -> SIC)
        [G_BS_users,I]=channelGain_BS(users_cordinates,alpha,sigma);
        G_BS_users=G_BS_users.^2;
        %generate BD-users channels before SIC order
        G_BD_users_unordered=channelGain_BD(BD_cordinates,users_cordinates,alpha);
        %ordering channels BD-users
        G_BD_users=G_BD_users_unordered(I);
        
        %% calling the function that computes R according to equation (6) in the paper
        R = rho_plus(sqrt(G_BS_users),G_BS_BD,G_BD_users);
        
        %% compute the minimum power budget (Pmin) required for meeting each QoS constraints in conventional NOMA and OMA
        for i=1:length(A)
            Pmin_OMA_conv(i)=sum((A(:,i).^K-1)./G_BS_users);
            
            Pmin_NOMA_conv(i)=(A(end,i)-1)/G_BS_users(end);
            for j=1:K
                Pmin_NOMA_conv(i)=Pmin_NOMA_conv(i)+(A(j,i)-1)/G_BS_users(j)*prod(A((j+1):K,i));
            end
            
            %% compute the optimal reflection coefficient(rho)
            if (isempty(R))
                rho_NOMA=1;
            else
                rho_NOMA=min(1,min(R));
            end
            
            %% compute Gamma according to the notations in the paper
            G_OMA_BD=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2;
            G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2;
            
            %% compute Pmin in OMA with BD and NOMA with BD
            Pmin_OMA_BD(i)=sum((A(i).^length(G_BS_users)-1)./G_OMA_BD);
            Pmin_NOMA_BD(i)=0;
            for j=1:length(G_BS_users)
                Pmin_NOMA_BD(i)=Pmin_NOMA_BD(i)+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:length(G_BS_users)));
            end
            
        end
    end
    
    %% stock results for each A (QoS constraints) and channel realization
    for i=1:length(A)
        %% compute the GEE for conventional NOMA/OMA and NOMA/OMA with backscatter
        EE_NOMA(n,i) = optimal_solution_NOMA(G_NOMA_BD,A(:,i),Pmax,Pmin_NOMA_BD(i),Pc);
        EE_OMA(n,i)= optimal_solution_OMA(G_OMA_BD,A(:,i),Pmax,Pc);
        EE_NOMA_conv(n,i) = optimal_solution_NOMA(G_BS_users,A(:,i),Pmax,Pmin_NOMA_conv(i),Pc);
        EE_OMA_conv(n,i) = optimal_solution_OMA(G_BS_users,A(:,i),Pmax,Pc);
        
    end
end

%% plot the figure (by averaging results over all channel realizations)

plot(Rmin,mean(EE_NOMA),'-o','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(Rmin));
hold on;
plot(Rmin,mean(EE_NOMA_conv),'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(Rmin));
hold on;
plot(Rmin,mean(EE_OMA),'-s','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(Rmin));
hold on;
plot(Rmin,mean(EE_OMA_conv),'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(Rmin));
ylabel('GEE (bits/J)');
xlabel('R_{min} (bits/s)');
legend('NOMA+backscatter','NOMA','OMA+backscatter','OMA','Location=Best');
grid on;
