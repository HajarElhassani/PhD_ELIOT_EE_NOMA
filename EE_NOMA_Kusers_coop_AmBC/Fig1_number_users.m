clear all;
close all;
clc;


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


%% compute the GEE as a function of the number of users K

K=2:10;%number of users

%for each channel realization
for n=1:N
    %for each user
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
        R = rho_plus(sqrt(G_BS_users),G_BS_BD,G_BD_users);
        
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
            R = rho_plus(sqrt(G_BS_users),G_BS_BD,G_BD_users);
            
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
        
        %% compute the optimal energy efficiency
        EE_NOMA_BD = optimal_solution_NOMA(G_NOMA_BD,A,Pmax,Pmin_NOMA_BD,Pc);
        EE_NOMA_conv = optimal_solution_NOMA(G_BS_users,A,Pmax,Pmin_NOMA_conv,Pc);
        EE_OMA_BD = optimal_solution_OMA(G_OMA_BD,A,Pmax,Pc);
        EE_OMA_conv = optimal_solution_OMA(G_BS_users,A,Pmax,Pc);
        
        
        %% stock results for each user
        EE_op_NOMA_BD(i)=EE_NOMA_BD;
        EE_op_NOMA_conv(i)=EE_NOMA_conv;
        EE_op_OMA_BD(i)=EE_OMA_BD;
        EE_op_OMA_conv(i)=EE_OMA_conv;
        
    end
    %% stock results for each channel realization 
    EE_op_NOMA_BD_n(n,:)=EE_op_NOMA_BD;
    EE_op_NOMA_conv_n(n,:)=EE_op_NOMA_conv;
    EE_op_OMA_BD_n(n,:)=EE_op_OMA_BD;
    EE_op_OMA_conv_n(n,:)=EE_op_OMA_conv;
    
end

%% averaging over channel realizations 
EE_op_NOMA_BD_mean=mean(EE_op_NOMA_BD_n);
EE_op_NOMA_conv_mean=mean(EE_op_NOMA_conv_n);
EE_op_OMA_BD_mean=mean(EE_op_OMA_BD_n);
EE_op_OMA_conv_mean=mean(EE_op_OMA_conv_n);

%% plot figures
figure(1)
plot(K,EE_op_NOMA_BD_mean,'-o','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_NOMA_conv_mean,'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_OMA_BD_mean,'-s','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_OMA_conv_mean,'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:length(K));
ylabel('GEE (bits/J)');
xlabel('Number of receivers');
legend('NOMA+backscatter','NOMA','OMA+backscatter','OMA','Location=Best');
grid on;
