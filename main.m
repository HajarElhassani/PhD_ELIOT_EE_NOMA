clear all;
clc;

%%Input parameters 
N=10^4; % number of channel realizations
channel_var=5; % users channel variance
Rmin=2; % QoS requirement
A=4*ones(1,3); % A=2^(Rmin)
Pmax_alpha=10; %the power budget at the BS in Watt
Pc=1; %the circuit power in Watt

%different values for the tradeoff parameters \alpha
alpha1=0; 
alpha2=0.2;
alpha3=0.5;
alpha4=1;
alpha5=1.5;

%% power allocation algorithm (NOMA/OMA)

%generate channel realizations
for k=1:N
    %for each channel realization
    G=channel_gain(3,channel_var); %generate the channel gains for the 3 users
    Pmin_OMA=((A(1)^3-1)/G(1)+(A(2)^3-1)/G(2)+(A(3)^3-1)/G(3))/3; %compute the minimum power required for OMA to satisfy QoS
    Pmin_NOMA=A(2)*A(3)*(A(1)-1)/G(1)+A(3)*(A(2)-1)/G(2)+(A(3)-1)/G(3); %compute the minimum power required for NOMA to satisfy QoS
    Pmin=max(Pmin_OMA,Pmin_NOMA); %the minimum power requirement for both OMA and NOMA
    % checking the feasibility condition
    while (Pmin_NOMA>Pmax_alpha||(A(1)^3-1)/G(1)>Pmax_alpha||(A(2)^3-1)/G(2)>Pmax_alpha||(A(3)^3-1)/G(3)>Pmax_alpha)
        G=channel_gain(3,channel_var);
        Pmin_OMA=((A(1)^3-1)/G(1)+(A(2)^3-1)/G(2)+(A(3)^3-1)/G(3))/3;
        Pmin_NOMA=A(2)*A(3)*(A(1)-1)/G(1)+A(3)*(A(2)-1)/G(2)+(A(3)-1)/G(3);
        Pmin=max(Pmin_OMA,Pmin_NOMA);
    end
    
    %%%%%%% part corresponding to figures 2a and 2b %%%%%%%%
    
    %%compute the relative gain and power excess w.r.t Pmax for each value of \alpha
    Pmax=(Pmin:1:Pmin+20); % the range of the power
    for i=1:length(Pmax) 
        [R1,PC1] = optimal_power_NOMA(A,G,alpha1,Pmax(i),Pc); %compute the optimal power allocation for NOMA
        [R_O1,PC_O1] = optimal_power_OMA(A,G,alpha1,Pmax(i),Pc); %compute the optimal power allocation for OMA
        SG1(i)=100*(R1-R_O1)/R_O1; % compute the relative gain
        PG1(i)=100*(PC_O1-PC1)/PC1; %compute the power excess
     
        [R2,PC2] = optimal_power_NOMA(A,G,alpha2,Pmax(i),Pc);
        [R_O2,PC_O2] = optimal_power_OMA(A,G,alpha2,Pmax(i),Pc);
        SG2(i)=100*(R2-R_O2)/R_O2;
        PG2(i)=100*(PC_O2-PC2)/PC2;
        
        [R3,PC3] = optimal_power_NOMA(A,G,alpha3,Pmax(i),Pc);
        [R_O3,PC_O3] = optimal_power_OMA(A,G,alpha3,Pmax(i),Pc);
        SG3(i)=100*(R3-R_O3)/R_O3;
        PG3(i)=100*(PC_O3-PC3)/PC3;
        
        [R4,PC4] = optimal_power_NOMA(A,G,alpha4,Pmax(i),Pc);
        [R_O4,PC_O4] = optimal_power_OMA(A,G,alpha4,Pmax(i),Pc);
        SG4(i)=100*(R4-R_O4)/R_O4;
        PG4(i)=100*(PC_O4-PC4)/PC4;
        
        [R5,PC5] = optimal_power_NOMA(A,G,alpha5,Pmax(i),Pc);
        [R_O5,PC_O5] = optimal_power_OMA(A,G,alpha5,Pmax(i),Pc);
        SG5(i)=100*(R5-R_O5)/R_O5;
        PG5(i)=100*(PC_O5-PC5)/PC5;
        
        Power_BS(i)=Pmax(i);   
    end
    
    %store each alpha(alpha1,...,alpha5) results vector above for each channel realization k
    SG_k1(k,:)=SG1;
    PG_k1(k,:)=PG1;
  
    SG_k2(k,:)=SG2;
    PG_k2(k,:)=PG2;
    
    SG_k3(k,:)=SG3;
    PG_k3(k,:)=PG3;
   
    SG_k4(k,:)=SG4;
    PG_k4(k,:)=PG4;
    
    SG_k5(k,:)=SG5;
    PG_k5(k,:)=PG5;
    
    Pmax_k(k,:)=Power_BS;
    
    %%%%%%% part corresponding to figures 1a and 1b %%%%%%%%
    
    %%compute the relative gain and power excess w.r.t alpha for fixed
    %%Pmax_alpha
    alpha=(0:.01:1.5); %the range of alpha
    for i=1:length(alpha)        
        [R,PC] = optimal_power_NOMA(A,G,alpha(i),Pmax_alpha,Pc);
        [R_O,PC_O] = optimal_power_OMA(A,G,alpha(i),Pmax_alpha,Pc);
        RG(i)=100*(R-R_O)/R_O;
        PG(i)=100*(PC_O-PC)/PC;
        
        %SG1_alpha(i)=RG;
        %PG1_alpha(i)=PG;
        
        %alpha_i(i)=alpha;
        %i=i+1;
        
    end
    
    SG_k(k,:)=RG;
    PG_k(k,:)=PG;
   
end

%% take the average of results over all channel realizations 

% results for figures 2a and 2b
SG1_mean=mean(SG_k1);
PG1_mean=mean(PG_k1);

SG2_mean=mean(SG_k2);
PG2_mean=mean(PG_k2);

SG3_mean=mean(SG_k3);
PG3_mean=mean(PG_k3);

SG4_mean=mean(SG_k4);
PG4_mean=mean(PG_k4);

SG5_mean=mean(SG_k5);
PG5_mean=mean(PG_k5);

Pmax_tot=round(mean(Pmax_k));


%results for figures 1a and 1b
SG_alpha_mean=mean(SG_k);
PG_alpha_mean=mean(PG_k);

%% plot figures

%figure 1a
figure(1);
plot(alpha,SG_alpha_mean,'LineWidth',2);
grid on;
xlabel('\alpha');
ylabel('G (%)');
%ffigure 1b
figure(2);
plot(alpha,PG_alpha_mean,'LineWidth',2);
grid on;
xlabel('\alpha');
ylabel('E (%)');

%figure 2a
figure(3);
plot(Pmax_tot,SG1_mean,'-o','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
hold on;
plot(Pmax_tot,SG2_mean,'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
hold on;
plot(Pmax_tot,SG3_mean,'-s','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
hold on;
plot(Pmax_tot,SG4_mean,'->','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
hold on;
plot(Pmax_tot,SG5_mean,'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
grid on;
xlabel('P_{max}');
ylabel('G(%)')
hLegend=legend({'\alpha=0','\alpha=0.2','\alpha=0.5','\alpha=1','\alpha=1.5'},'Location','Best');
set(hLegend, 'Color','none');

%figure 2b
figure(4);
plot(Pmax_tot,PG1_mean,'-o','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
hold on;
plot(Pmax_tot,PG2_mean,'-*','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
hold on;
plot(Pmax_tot,PG3_mean,'-s','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
hold on;
plot(Pmax_tot,PG4_mean,'->','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
hold on;
plot(Pmax_tot,PG5_mean,'-d','MarkerSize',7,'LineWidth',1.3,'MarkerIndices', 1:2:length(Pmax_tot));
grid on;
xlabel('P_{max}');
ylabel('E (%)');
hLegend=legend({'\alpha=0','\alpha=0.2','\alpha=0.5','\alpha=1','\alpha=1.5'},'Location','Best');
set(hLegend, 'Color','none');


