clear all;
close all;
clc;


%---------- coordinates to create a cell ----------%


%the minimum x and y coordinates of BD (backscatter device) from BS
b1=0;
%the maximum x and y coordinates of BD from BS
b2=3;

%the minimum x and y coordinates of users from BS
a1=0;
%the maximum x and y coordinates of users from BS
a2=15;

%---------- input parameters ----------%


%% Input parameters

N=1e3;%number of channel realizations
Pmax_dbm=80;%power budget at the BS in dbm
Pc_dbm=30;%circuit power in dbm
sigma_dbm=-20;%noise power in dbm
Rmin=1;%minimum  rate QoS constraint

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

%% compute rho as a function of number of users

K=2:10;%number of users


%% for each channel realization
for n=1:N
    %% for each number of users K
    for i=1:length(K)
        A=(2^(2*Rmin))*ones(K(i),1);
        
        %% generating channels
        %generate x and y coordinates for users
        users_cordinates = coordinates(K(i),radius_BS_users,min_dis_BS_users)';
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
        R = rho_plus(G_BS_users,G_BS_BD,G_BD_users);
        
        %% compute the minimum power budget (Pmin) required for meeting QoS constraints in conventional NOMA
        Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K(i)));
        end
        
        %% compute Pmin required for meeting QoS constraints in conventional OMA
        Pmin_OMA_conv=sum((A.^K(i)-1)./G_BS_users);
        
        %% compute the optimal rho
        if (isempty(R))
            rho_NOMA=1;
        else
            rho_NOMA=min(1,min(R));
        end
        
        %% compute the backscattered channel (Gamma) according to the notations in the paper
        G_OMA_BD=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2;
        G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2;
        
        %% compute Pmin in OMA with BD and NOMA with BD
        Pmin_OMA_BD=sum((A.^length(G_BS_users)-1)./G_OMA_BD);
        
        Pmin_NOMA_BD=(A(end)-1)/G_NOMA_BD(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD=Pmin_NOMA_BD+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:K(i)));
        end
        
        %% checking the feasability condition: if not satisfied repeat generating channel gains until the feasibility condition is met
        while (Pmin_NOMA_conv>Pmax || Pmin_NOMA_BD>Pmax || Pmin_OMA_BD>Pmax || Pmin_OMA_conv>Pmax)
            
            %generate x and y coordinates for users
            users_cordinates = coordinates(K(i),radius_BS_users,min_dis_BS_users)';
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
            R = rho_plus(G_BS_users,G_BS_BD,G_BD_users);
            
            %% compute the minimum power budget (Pmin) required for meeting QoS constraints in conventional NOMA
            Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K(i)));
            end
            
            %% compute Pmin required for meeting QoS constraints in conventional OMA
            Pmin_OMA_conv=sum((A.^K(i)-1)./G_BS_users);
            
            %% compute the optimal rho
            if (isempty(R))
                rho_NOMA=1;
            else
                rho_NOMA=min(1,min(R));
            end
            
            %% compute the backscattered channel (Gamma) according to the notations in the paper
            G_OMA_BD=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2;
            G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2;
            
            %% compute Pmin in OMA with BD and NOMA with BD
            Pmin_OMA_BD=sum((A.^length(G_BS_users)-1)./G_OMA_BD);
            
            Pmin_NOMA_BD=(A(end)-1)/G_NOMA_BD(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD=Pmin_NOMA_BD+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:K(i)));
            end
            
        end
        
        %% stock results for each number of users 
        rho_op(i)=rho_NOMA;
    end
    
    %% stock results for each channel realization 
    rho_op_n(n,:)=rho_op;
    
end

%% averaging over channel realizations
rho_op_mean=mean(rho_op_n);

%% plot figures
plot(K,rho_op_mean,'-o','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
hold on;
plot(K,1*ones(length(K),1),'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
ylabel('\rho^*');
xlabel('Number of receivers');
legend('NOMA+backscatter','OMA+backscatter','Location=Best');
grid on;