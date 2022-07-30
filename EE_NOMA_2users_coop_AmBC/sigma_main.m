clear all;
close all;
clc;

%%  Figure.2: Energy efficiency (GEE) as a function of sigma 

%% input parameters
min_dis_BS_users=0;%the minimum distance between BS and users
radius_BS_users=15;%radius of users cell

min_dis_BS_BD1=0;%the minimum distance between BD1 and BS
min_dis_BS_BD2=0;%the minimum distance between BD2 and BS
min_dis_BS_BD3=0;%the minimum distance between BD3 and BS
min_dis_BS_BD4=0;%the minimum distance between BD4 and BS

radius_BS_BD1=3;%radius of BD1 cell
radius_BS_BD2=3;%radius of BD2 cell
radius_BS_BD3=3;%radius of BD3 cell
radius_BS_BD4=3;%radius of BD4 cell

alpha_path=2.5; %pathloss exponent

sigma_dbm=-50:10:0;% the noise variance in dBm
Pmax_dbm=40; %power budget in dbm
Pc_dbm=30; %circuit power in dBm

%conversion from dBm to Watt
sigma=dbm_to_Watt(sigma_dbm);
Pmax=dbm_to_Watt(Pmax_dbm);
Pc=dbm_to_Watt(Pc_dbm);

K=2; %number of users

Rmin=1; %QoS constraints
A=(2^(2*Rmin))*ones(K,1);

N=10 ; %number of realizations

for n=1:N
    for i=1:length(sigma)
        
        %% generate channel gains
        %generate x and y coordinates for users within the cell
        users_cordinates = coordinates(K,radius_BS_users,min_dis_BS_users)';
        %generate x and y coordinates for all backscatter devices
        BD1_cordinates = coordinates(1,radius_BS_BD1,min_dis_BS_BD1);
        BD2_cordinates = coordinates(1,radius_BS_BD2,min_dis_BS_BD2);
        BD3_cordinates = coordinates(1,radius_BS_BD3,min_dis_BS_BD3);
        BD4_cordinates = coordinates(1,radius_BS_BD4,min_dis_BS_BD4);
        
        %generate BD-BS channels
        G_BS_BD1=channelGain_BS(BD1_cordinates,alpha_path,sigma(i));
        G_BS_BD2=channelGain_BS(BD2_cordinates,alpha_path,sigma(i));
        G_BS_BD3=channelGain_BS(BD3_cordinates,alpha_path,sigma(i));
        G_BS_BD4=channelGain_BS(BD4_cordinates,alpha_path,sigma(i));
        
        %generate BS-users channel (following SIC order 'descend')
        [G_users,I]=channelGain_BS(users_cordinates,alpha_path,sigma(i));
        
        %generate BD-users channels before SIC order
        G_BD1_users_unordered=channelGain_BD(BD1_cordinates,users_cordinates,alpha_path);
        G_BD2_users_unordered=channelGain_BD(BD2_cordinates,users_cordinates,alpha_path);
        G_BD3_users_unordered=channelGain_BD(BD3_cordinates,users_cordinates,alpha_path);
        G_BD4_users_unordered=channelGain_BD(BD4_cordinates,users_cordinates,alpha_path);
        
        %ordering BD-users channels
        G_BD1_users=G_BD1_users_unordered(I);
        G_BD2_users=G_BD2_users_unordered(I);
        G_BD3_users=G_BD3_users_unordered(I);
        G_BD4_users=G_BD4_users_unordered(I);
        
        %BS-users channel gain
        G_users=G_users.^2;
        
        %pick randomly one of the two backscatter devices to compare with NOMA with one backscatter device.
        idx = randperm(length([G_BS_BD1,G_BS_BD2]),1);
        if (idx==1)
            G_BS_BD=G_BS_BD1;
            G_BD_users=G_BD1_users;
        else
            G_BS_BD=G_BS_BD2;
            G_BD_users=G_BD2_users;
        end
        
        
        %% compute the optimal reflection coefficients and the channels of the backscattered links
        
        [rho_4BD,G_NOMA_4BD] = optimal_rho_4BD(sqrt(G_users),G_BS_BD1,G_BS_BD2,G_BS_BD3,G_BS_BD4,G_BD1_users,G_BD2_users,G_BD3_users,G_BD4_users);
        [rho_3BD,G_NOMA_3BD] = optimal_rho_3BD(sqrt(G_users),G_BS_BD1,G_BS_BD2,G_BS_BD3,G_BD1_users,G_BD2_users,G_BD3_users);
        [rho_2BD,G_NOMA_2BD,G_OMA_2BD] = optimal_rho_2BD(sqrt(G_users),G_BS_BD1,G_BS_BD2,G_BD1_users,G_BD2_users);
        [rho_1BD,G_NOMA_1BD] = optimal_rho_1BD(sqrt(G_users),G_BS_BD,G_BD_users);
        
        %% compute the minimum power for OMA and NOMA (comventional, backscattered links :1BD; 2BD; 3BD; 4BD)
        Pmin_NOMA_conv=(A(2)-1)/G_users(2)+A(2)*(A(1)-1)/G_users(1);
        Pmin_NOMA_1BD=(A(2)-1)/G_NOMA_1BD(2)+A(2)*(A(1)-1)/G_NOMA_1BD(1);
        Pmin_NOMA_2BD=(A(2)-1)/G_NOMA_2BD(2)+A(2)*(A(1)-1)/G_NOMA_2BD(1);
        Pmin_NOMA_3BD=(A(2)-1)/G_NOMA_3BD(2)+A(2)*(A(1)-1)/G_NOMA_3BD(1);
        Pmin_NOMA_4BD=(A(2)-1)/G_NOMA_4BD(2)+A(2)*(A(1)-1)/G_NOMA_4BD(1);
        %% minimum power for OMA (two backscatter devices) backscattered link
        Pmin_OMA_2BD=sum((A.^length(G_OMA_2BD)-1)./G_OMA_2BD);
        
        
        while (Pmin_NOMA_conv>Pmax || Pmin_NOMA_2BD>Pmax || Pmin_NOMA_1BD>Pmax || Pmin_OMA_2BD>Pmax)
            
            %% generate channel gains
            %generate x and y coordinates for users within the cell
            users_cordinates = coordinates(K,radius_BS_users,min_dis_BS_users)';
            %generate x and y coordinates for all backscatter devices
            BD1_cordinates = coordinates(1,radius_BS_BD1,min_dis_BS_BD1);
            BD2_cordinates = coordinates(1,radius_BS_BD2,min_dis_BS_BD2);
            BD3_cordinates = coordinates(1,radius_BS_BD3,min_dis_BS_BD3);
            BD4_cordinates = coordinates(1,radius_BS_BD4,min_dis_BS_BD4);
            
            %generate BD-BS channels
            G_BS_BD1=channelGain_BS(BD1_cordinates,alpha_path,sigma(i));
            G_BS_BD2=channelGain_BS(BD2_cordinates,alpha_path,sigma(i));
            G_BS_BD3=channelGain_BS(BD3_cordinates,alpha_path,sigma(i));
            G_BS_BD4=channelGain_BS(BD4_cordinates,alpha_path,sigma(i));
            
            %generate BS-users channel (following SIC order 'descend')
            [G_users,I]=channelGain_BS(users_cordinates,alpha_path,sigma(i));
            
            %generate BD-users channels before SIC order
            G_BD1_users_unordered=channelGain_BD(BD1_cordinates,users_cordinates,alpha_path);
            G_BD2_users_unordered=channelGain_BD(BD2_cordinates,users_cordinates,alpha_path);
            G_BD3_users_unordered=channelGain_BD(BD3_cordinates,users_cordinates,alpha_path);
            G_BD4_users_unordered=channelGain_BD(BD4_cordinates,users_cordinates,alpha_path);
            
            %ordering BD-users channels
            G_BD1_users=G_BD1_users_unordered(I);
            G_BD2_users=G_BD2_users_unordered(I);
            G_BD3_users=G_BD3_users_unordered(I);
            G_BD4_users=G_BD4_users_unordered(I);
            
            %BS-users channel gain
            G_users=G_users.^2;
            
            %pick randomly one of the two backscatter devices to compare with NOMA with one backscatter device.
            idx = randperm(length([G_BS_BD1,G_BS_BD2]),1);
            if (idx==1)
                G_BS_BD=G_BS_BD1;
                G_BD_users=G_BD1_users;
            else
                G_BS_BD=G_BS_BD2;
                G_BD_users=G_BD2_users;
            end
            
            
            %% compute the optimal reflection coefficients and the channels of the backscattered links
            
            [rho_4BD,G_NOMA_4BD] = optimal_rho_4BD(sqrt(G_users),G_BS_BD1,G_BS_BD2,G_BS_BD3,G_BS_BD4,G_BD1_users,G_BD2_users,G_BD3_users,G_BD4_users);
            [rho_3BD,G_NOMA_3BD] = optimal_rho_3BD(sqrt(G_users),G_BS_BD1,G_BS_BD2,G_BS_BD3,G_BD1_users,G_BD2_users,G_BD3_users);
            [rho_2BD,G_NOMA_2BD,G_OMA_2BD] = optimal_rho_2BD(G_users,G_BS_BD1,G_BS_BD2,G_BD1_users,G_BD2_users);
            [rho_1BD,G_NOMA_1BD] = optimal_rho_1BD(G_users,G_BS_BD,G_BD_users);
            
            %% compute the minimum power for OMA and NOMA (comventional, backscattered links :1BD; 2BD; 3BD; 4BD)
            Pmin_NOMA_conv=(A(2)-1)/G_users(2)+A(2)*(A(1)-1)/G_users(1);
            Pmin_NOMA_1BD=(A(2)-1)/G_NOMA_1BD(2)+A(2)*(A(1)-1)/G_NOMA_1BD(1);
            Pmin_NOMA_2BD=(A(2)-1)/G_NOMA_2BD(2)+A(2)*(A(1)-1)/G_NOMA_2BD(1);
            Pmin_NOMA_3BD=(A(2)-1)/G_NOMA_3BD(2)+A(2)*(A(1)-1)/G_NOMA_3BD(1);
            Pmin_NOMA_4BD=(A(2)-1)/G_NOMA_4BD(2)+A(2)*(A(1)-1)/G_NOMA_4BD(1);
            %% minimum power for OMA (two backscatter devices) backscattered link
            Pmin_OMA_2BD=sum((A.^length(G_OMA_2BD)-1)./G_OMA_2BD);
            
        end
        
        %% compute NOMA energy efficiency ratio for each sigma and channel realization
        EE_NOMA_2BD(i) = optimal_solution(G_NOMA_2BD,A,Pmax,Pmin_NOMA_2BD,Pc);
        EE_NOMA_1BD(i) = optimal_solution_NOMA(G_NOMA_1BD,A,Pmax,Pmin_NOMA_1BD,Pc);
        EE_NOMA_conv(i) = optimal_solution_NOMA(G_users,A,Pmax,Pmin_NOMA_conv,Pc);
        EE_OMA_2BD(i) = optimal_solution_OMA(G_OMA_2BD,A,Pmax,Pc);
        EE_NOMA_3BD(i) = optimal_solution_NOMA(G_NOMA_3BD,A,Pmax,Pmin_NOMA_3BD,Pc);
        EE_NOMA_4BD(i) = optimal_solution_NOMA(G_NOMA_4BD,A,Pmax,Pmin_NOMA_4BD,Pc);
        
    end
    
    %% compute the relative gain w.r.t conventional NOMA
    R_EE_NOMA_1BD(n,:)=100*(EE_NOMA_1BD-EE_NOMA_conv)./EE_NOMA_conv;
    R_EE_NOMA_2BD(n,:)=100*(EE_NOMA_2BD-EE_NOMA_conv)./EE_NOMA_conv;
    R_EE_NOMA_3BD(n,:)=100*(EE_NOMA_3BD-EE_NOMA_conv)./EE_NOMA_conv;
    R_EE_NOMA_4BD(n,:)=100*(EE_NOMA_4BD-EE_NOMA_conv)./EE_NOMA_conv;
    R_EE_OMA_2BD(n,:)=100*(EE_OMA_2BD-EE_NOMA_conv)./EE_NOMA_conv; 
    
    if (mod(n,100)==0)
        disp(n)
    end
    
    
end

%% averaging over all the channel realizations
R_EE_NOMA_2BD_mean=mean(R_EE_NOMA_2BD);
R_EE_OMA_2BD_mean=mean(R_EE_OMA_2BD);
R_EE_NOMA_1BD_mean=mean(R_EE_NOMA_1BD);
R_EE_NOMA_3BD_mean=mean(R_EE_NOMA_3BD);
R_EE_NOMA_4BD_mean=mean(R_EE_NOMA_4BD);

%% plot figure
figure(1)
plot(sigma_dbm,mean(R_EE_NOMA_2BD),'-o','MarkerSize',7,'LineWidth',2.5);
hold on;
plot(sigma_dbm,mean(R_EE_OMA_2BD),'-s','MarkerSize',7,'LineWidth',2.5);
hold on;
plot(sigma_dbm,mean(R_EE_NOMA_1BD),'-*','MarkerSize',7,'LineWidth',2.5);
hold on;
plot(sigma_dbm,mean(R_EE_NOMA_3BD),':|','MarkerSize',7,'LineWidth',2.5);
hold on;
plot(sigma_dbm,mean(R_EE_NOMA_4BD),':|','MarkerSize',7,'LineWidth',2.5);
ylabel('Relative gain (%)');
xlabel('\sigma^2 (dBm)');
legend('NOMA+2 backscatters vs NOMA','OMA+2 backscatters vs NOMA','NOMA+1 backscatter vs NOMA','NOMA+3 backscatters vs NOMA','NOMA+4 backscatters vs NOMA','Location=Best');
grid on;
