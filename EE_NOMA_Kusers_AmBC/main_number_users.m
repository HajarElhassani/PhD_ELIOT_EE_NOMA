clear all;
close all;
clc;


%% Input parameters
N=1e3;%number of realizations
Rmin=1;%minimum  rate QoS constraint

Pmax_dbm=60;%power budget in dbm
Pc_dbm=30; %circuit power in dbm
sigma_dbm=-20;%noise power in dbm

%call the function 'dbm_to_Watt' to convert from dbm to Watt
Pmax=dbm_to_Watt(Pmax_dbm);
Pc=dbm_to_Watt(Pc_dbm);
sigma=dbm_to_Watt(sigma_dbm);

%probability of backscattering (B=1)
q1=1;
q2=1/2;
q3=0;

%create the cell
min_dis_BS_users=0;%the minimum distance between BS and users
radius_BS_users=20;%maximum distance between BS and users
min_dis_BS_BD=0;%the minimum distance between BD and BS
radius_BS_BD=4;%maximum distance between BD and BS
alpha=2.5;%pathloss exponent

%-------- compute the GEE as a function of the number of users K ----------%
K=2:7;%number of users

%for each channel realization
for n=1:N
    %for each number of users K
    for i=1:length(K)
        % SNR threshold
        A=(2^(2*Rmin))*ones(K(i),1);
        
        %---- generating channel gains ----%
        %generate x and y coordinates for users
        users_cordinates = coordinates(K(i),radius_BS_users,min_dis_BS_users)';
        %generate x and y coordinates for BD
        BD_cordinates = coordinates(1,radius_BS_BD,min_dis_BS_BD);
        
        %generate BS-BD channel gain
        G_BS_BD=channelGain_BS(1,BD_cordinates,alpha,sigma,0);
        %generate BS-users channel gains (in descending order -> SIC)
        [G_BS_users,I]=channelGain_BS(1,users_cordinates,alpha,sigma,0);
        G_BS_users=G_BS_users.^2; % channel gain
        %generate BD-users channel gains before SIC order
        G_BD_users_unordered=channelGain_BD(BD_cordinates,users_cordinates,alpha,0);
        %ordering channel gains BD-users
        G_BD_users=G_BD_users_unordered(I);
        
        
        
        %---- calling the function that computes R according to equation (6) in the paper ----%
        R = rho_plus(sqrt(G_BS_users),G_BS_BD,G_BD_users);
        
        %---- compute the minimum power budget (Pmin) required for meeting QoS constraints in conventional NOMA ----%
        Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K(i)));
        end
        
        %---- compute Pmin required for meeting QoS constraints in conventional OMA ----%
        Pmin_OMA_conv=sum((A.^K(i)-1)./G_BS_users);
        
        %compute the optimal reflection coefficient rho
        if (isempty(R))
            rho_NOMA=1;
        else
            rho_NOMA=min(1,min(R));
        end
        
        %compute (the backscattered+direct channel gain) Gamma according to the notations in the paper
        G_OMA_BD=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2;
        G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2;
        
        %compute Pmin in OMA with BD and NOMA with BD
        Pmin_OMA_BD=sum((A.^length(G_BS_users)-1)./G_OMA_BD);
        
        Pmin_NOMA_BD=(A(end)-1)/G_NOMA_BD(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD=Pmin_NOMA_BD+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:K(i)));
        end
        
        
        %% checking the feasability condition
        
        while (Pmin_NOMA_conv>Pmax || Pmin_NOMA_BD>Pmax || Pmin_OMA_BD>Pmax || Pmin_OMA_conv>Pmax)
            
            %generate x and y coordinates for users
            users_cordinates = coordinates(K(i),radius_BS_users,min_dis_BS_users)';
            %generate x and y coordinates for BD
            BD_cordinates = coordinates(1,radius_BS_BD,min_dis_BS_BD);
            
            %generate BS-BD channel gain
            G_BS_BD=channelGain_BS(1,BD_cordinates,alpha,sigma,0);
            %generate BS-users channel gains (in descending order -> SIC)
            [G_BS_users,I]=channelGain_BS(1,users_cordinates,alpha,sigma,0);
            G_BS_users=G_BS_users.^2; % channel gain
            %generate BD-users channel gains before SIC order
            G_BD_users_unordered=channelGain_BD(BD_cordinates,users_cordinates,alpha,0);
            %ordering channel gains BD-users
            G_BD_users=G_BD_users_unordered(I);
            
            
            
            %---- calling the function that computes R according to equation (6) in the paper ----%
            R = rho_plus(sqrt(G_BS_users),G_BS_BD,G_BD_users);
            
            %---- compute the minimum power budget (Pmin) required for meeting QoS constraints in conventional NOMA ----%
            Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K(i)));
            end
            
            %---- compute Pmin required for meeting QoS constraints in conventional OMA ----%
            Pmin_OMA_conv=sum((A.^K(i)-1)./G_BS_users);
            
            %compute the optimal reflection coefficient rho
            if (isempty(R))
                rho_NOMA=1;
            else
                rho_NOMA=min(1,min(R));
            end
            
            %compute (the backscattered+direct channel gain) Gamma according to the notations in the paper
            G_OMA_BD=(sqrt(G_BS_users)+G_BS_BD*G_BD_users).^2;
            G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2;
            
            %compute Pmin in OMA with BD and NOMA with BD
            Pmin_OMA_BD=sum((A.^length(G_BS_users)-1)./G_OMA_BD);
            
            Pmin_NOMA_BD=(A(end)-1)/G_NOMA_BD(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD=Pmin_NOMA_BD+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:K(i)));
            end
            
        end
        
        %---- compute the optimal energy efficiency for different schemes and stock results----%
        
        EE_NOMA_BD1(i) = Dinkelbach_NOMA(q1,G_BS_users,G_NOMA_BD,A,Pmax,Pmin_NOMA_conv,Pc);
        EE_NOMA_BD2(i) = Dinkelbach_NOMA(q2,G_BS_users,G_NOMA_BD,A,Pmax,Pmin_NOMA_conv,Pc);
        EE_NOMA_BD3(i) = optimal_solution_NOMA(G_BS_users,A,Pmax,Pmin_NOMA_conv,Pc);
        EE_OMA_BD1(i) = Dinkelbach_OMA(q1,G_BS_users,G_OMA_BD,A,Pmax,Pc);
        EE_OMA_BD2(i) = Dinkelbach_OMA(q2,G_BS_users,G_OMA_BD,A,Pmax,Pc);
        EE_OMA_BD3(i) = optimal_solution_OMA(G_BS_users,A,Pmax,Pc);        
        
    end
    %-------- stock results for each channel realization --------%
    EE_op_NOMA_BD_n1(n,:)=EE_NOMA_BD1;
    EE_op_NOMA_BD_n2(n,:)=EE_NOMA_BD2;
    EE_op_NOMA_BD_n3(n,:)=EE_NOMA_BD3;
    EE_op_OMA_BD_n1(n,:)=EE_OMA_BD1;
    EE_op_OMA_BD_n2(n,:)=EE_OMA_BD2;
    EE_op_OMA_BD_n3(n,:)=EE_OMA_BD3;
    
    % the relative gain
    R_NOMA_OMA_BD1(n,:)=100*(EE_NOMA_BD1-EE_OMA_BD1)./EE_OMA_BD1;
    R_NOMA_OMA_BD2(n,:)=100*(EE_NOMA_BD2-EE_OMA_BD1)./EE_OMA_BD1;
    
    R_NOMA_convNOMA_BD1(n,:)=100*(EE_NOMA_BD1-EE_NOMA_BD3)./EE_NOMA_BD3;
    R_NOMA_convNOMA_BD2(n,:)=100*(EE_NOMA_BD2-EE_NOMA_BD3)./EE_NOMA_BD3;
    
    
end

%------------- averaging over channel realizations -----------%

%energy efficiency
EE_op_NOMA_BD_mean1=mean(EE_op_NOMA_BD_n1);
EE_op_NOMA_BD_mean2=mean(EE_op_NOMA_BD_n2);
EE_op_NOMA_BD_mean3=mean(EE_op_NOMA_BD_n3);

EE_op_OMA_BD_mean1=mean(EE_op_OMA_BD_n1);
EE_op_OMA_BD_mean2=mean(EE_op_OMA_BD_n2);
EE_op_OMA_BD_mean3=mean(EE_op_OMA_BD_n3);

%relative gains
R_NOMA_OMA_BD1_mean=mean(R_NOMA_OMA_BD1);
R_NOMA_OMA_BD2_mean=mean(R_NOMA_OMA_BD2);

R_NOMA_convNOMA_BD1_mean=mean(R_NOMA_convNOMA_BD1);
R_NOMA_convNOMA_BD2_mean=mean(R_NOMA_convNOMA_BD2);

%------------- plot figures -----------%
figure(1)
plot(K,EE_op_NOMA_BD_mean1,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_NOMA_BD_mean2,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_NOMA_BD_mean3,'-*','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_OMA_BD_mean1,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_OMA_BD_mean2,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_op_OMA_BD_mean3,'-*','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
ylabel('\xi_{EE} (bits/J)');
xlabel('Number of receivers');
legend('NOMA(q=1)','NOMA(q=0.5)','NOMA(q=0)','OMA (q=1)','OMA (q=0.5)','OMA (q=0)');
grid on;

figure(2)
plot(K,R_NOMA_OMA_BD1_mean,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,R_NOMA_OMA_BD2_mean,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
ylabel('relative gain to OMA (q=1)');
xlabel('Number of receivers');
legend('NOMA(q=1)','NOMA(q=0.5)');
grid on;

figure(3)
plot(K,R_NOMA_convNOMA_BD2_mean,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
ylabel('relative gain to conv NOMA');
xlabel('Number of receivers');
legend('NOMA(q=0.5)');
grid on;