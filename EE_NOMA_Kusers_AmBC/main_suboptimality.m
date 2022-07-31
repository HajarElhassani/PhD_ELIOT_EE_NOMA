clear all;
close all;
clc;

%% Input parameters
N=1e3;%number of
K=2; %number of receivers
Rmin=1; %minimum  rate QoS constraint
A=(2^(2*Rmin))*ones(K,1);

Pmax_dbm=30; %power budget in dbm
Pc_dbm=30; %circuit power in dbm
sigma_dbm=-20; %noise power in dbm

%call the function 'dbm_to_Watt' to convert from dbm to Watt
Pmax=dbm_to_Watt(Pmax_dbm);
Pc=dbm_to_Watt(Pc_dbm);
sigma=dbm_to_Watt(sigma_dbm);


%coordinates to create the cell
min_dis_BS_users=0;%the minimum distance between BS and users
radius_BS_users=20;%maximum distance between BS and users
min_dis_BS_BD=0;%the minimum distance between BD and BS
radius_BS_BD=4;%maximum distance between BD and BS
alpha=2.5; %pathloss exponent

q=0:.2:1; %probability of backscattering state (B=1)

%% compare our solution with the optimal one obtained by exhaustive search for K=2

%for each channel realization
for n=1:N
    
    %% generating channel gains
    %generate x and y coordinates for users
    users_cordinates = coordinates(K,radius_BS_users,min_dis_BS_users)';
    %generate x and y coordinates for BD
    BD_cordinates = coordinates(1,radius_BS_BD,min_dis_BS_BD);
    
    %generate BS-BD channel gain
    G_BS_BD=channelGain_BS(1,BD_cordinates,alpha,sigma,0);
    %generate BS-users channel gains (in descending order -> SIC)
    [G_BS_users,I]=channelGain_BS(1,users_cordinates,alpha,sigma,0);
    G_BS_users=G_BS_users.^2;
    %generate BD-users channel gains before SIC order
    G_BD_users_unordered=channelGain_BD(BD_cordinates,users_cordinates,alpha,0);
    %ordering channel gains BD-users
    G_BD_users=G_BD_users_unordered(I);
    
    %---- calling the function that computes R according to equation (6) in the paper ----%
    R = rho_plus(sqrt(G_BS_users),G_BS_BD,G_BD_users);
    
    %---- compute the minimum power budget (Pmin) required for meeting QoS constraints in conventional NOMA ----%
    Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
    for j=1:(K-1)
        Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K));
    end
    
    %---- compute Pmin required for meeting QoS constraints in conventional OMA ----%
    Pmin_OMA_conv=sum((A.^K-1)./G_BS_users);
    
    %compute the optimal rho
    if (isempty(R))
        rho_NOMA=1;
    else
        rho_NOMA=min(1,min(R));
    end
    
    %compute the backscattered+direct channel gains (Gamma) according to the notations in the paper
    G_OMA_BD=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2; % for OMA
    G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2; % for NOMA
    
    %compute Pmin in OMA with BD and NOMA with BD
    Pmin_OMA_BD=sum((A.^length(G_BS_users)-1)./G_OMA_BD);
    
    Pmin_NOMA_BD=(A(end)-1)/G_NOMA_BD(end);
    for j=1:(K-1)
        Pmin_NOMA_BD=Pmin_NOMA_BD+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:K));
    end
    
    %% checking the feasability condition
    while (Pmin_NOMA_conv>Pmax || Pmin_NOMA_BD>Pmax || Pmin_OMA_BD>Pmax || Pmin_OMA_conv>Pmax)
        
        %generate x and y coordinates for users
        users_cordinates = coordinates(K,radius_BS_users,min_dis_BS_users)';
        %generate x and y coordinates for BD
        BD_cordinates = coordinates(1,radius_BS_BD,min_dis_BS_BD);
        
        %generate BS-BD channel gain
        G_BS_BD=channelGain_BS(1,BD_cordinates,alpha,sigma,0);
        %generate BS-users channel gains (in descending order -> SIC)
        [G_BS_users,I]=channelGain_BS(1,users_cordinates,alpha,sigma,0);
        G_BS_users=G_BS_users.^2;
        %generate BD-users channel gains before SIC order
        G_BD_users_unordered=channelGain_BD(BD_cordinates,users_cordinates,alpha,0);
        %ordering channel gains BD-users
        G_BD_users=G_BD_users_unordered(I);
        
        %---- calling the function that computes R according to equation (6) in the paper ----%
        R = rho_plus(sqrt(G_BS_users),G_BS_BD,G_BD_users);
        
        %---- compute the minimum power budget (Pmin) required for meeting QoS constraints in conventional NOMA ----%
        Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
        for j=1:(K-1)
            Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K));
        end
        
        %---- compute Pmin required for meeting QoS constraints in conventional OMA ----%
        Pmin_OMA_conv=sum((A.^K-1)./G_BS_users);
        
        %compute the optimal rho
        if (isempty(R))
            rho_NOMA=1;
        else
            rho_NOMA=min(1,min(R));
        end
        
        %compute the backscattered+direct channel gains (Gamma) according to the notations in the paper
        G_OMA_BD=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2; % for OMA
        G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2; % for NOMA
        
        %compute Pmin in OMA with BD and NOMA with BD
        Pmin_OMA_BD=sum((A.^length(G_BS_users)-1)./G_OMA_BD);
        
        Pmin_NOMA_BD=(A(end)-1)/G_NOMA_BD(end);
        for j=1:(K-1)
            Pmin_NOMA_BD=Pmin_NOMA_BD+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:K));
        end
        
    end
    
    %---- compute the optimal energy efficiency for every probability q ----%
    
    %for each probability
    for i=1:length(q) 
        EE_NOMA_sub(i) = Dinkelbach_NOMA(q(i),G_BS_users,G_NOMA_BD,A,Pmax,Pmin_NOMA_conv,Pc);
        EE_NOMA_op(i) = exhaustive_solution(q(i),sqrt(G_BS_users),G_BS_BD*G_BD_users,Pmax,Pc,Rmin);
    end
    %-------- stock results for each channel realization --------%
    EE_NOMA_sub_n(n,:)=EE_NOMA_sub;
    EE_NOMA_op_n(n,:)=EE_NOMA_op;
    
end

%------------- averaging over channel realizations -----------%
EE_NOMA_sub_mean=mean(EE_NOMA_sub_n);
EE_NOMA_op_mean=mean(EE_NOMA_op_n);

%------------- plot figures -----------%
figure(1)
plot(q,EE_NOMA_sub_mean,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(q));
hold on;
plot(q,EE_NOMA_op_mean,'-*','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(q));
ylabel('\xi_{EE} (bits/J)');
xlabel('q');
legend('NOMA sub- optimal', 'NOMA optimal');
grid on;
