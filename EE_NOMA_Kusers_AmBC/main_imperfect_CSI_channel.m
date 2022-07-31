clear all;
close all;
clc;

%% input parameters

N=1e2;%number of realizations
Rmin=1;%minimum  rate QoS constraint
q=1/2;% probability of backscattering (B=1)

Pmax_dbm=60;%power budget in dbm
Pc_dbm=30;%circuit power in dbm
sigma_dbm=-20;%noise power in dbm

%call the function 'dbm_to_Watt' to convert from dbm to Watt
Pmax=dbm_to_Watt(Pmax_dbm);
Pc=dbm_to_Watt(Pc_dbm);
sigma=dbm_to_Watt(sigma_dbm);


%create the cell
min_dis_BS_users=0;%the minimum distance between BS and users
radius_BS_users=20;%maximum distance between BS and users
min_dis_BS_BD=0;%the minimum distance between BD and BS
radius_BS_BD=4;%maximum distance between BD and BS
alpha=2.5;%pathloss exponent


%error variance computed via montecarlo such that SNR = {20,10,0,-10,-20}
%dB
error1=[0.000747,1e-9,0]; %20dB
error2=[0.00747,1e-8,0]; %10dB
error3=[0.0747,1e-7,0]; %0dB
error4=[0.747,1e-6,0]; %-10dB
error5=[7.47,1e-5,0]; %-20dB

%-------- compute the GEE in case of imperfect channel as a function of the number of users K ----------%
K=2:7;%number of users

%for each channel realization
for n=1:N
    for i=1:length(K)
        % the SNR threshold
        A=(2^(2*Rmin))*ones(K(i),1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE OF PERFECT CSI %%%%%%%%%%%%%%%
        %generate x and y coordinates for users
        users_cordinates = coordinates(K(i),radius_BS_users,min_dis_BS_users)';
        %generate x and y coordinates for BD
        BD_cordinates = coordinates(1,radius_BS_BD,min_dis_BS_BD);
        
        %generate BS-BD channel gain
        G_BS_BD=channelGain_BS(0,BD_cordinates,alpha,sigma,0);
        %generate BS-users channel gains (in descending order -> SIC)
        [G_BS_users_unordered,I]=channelGain_BS(0,users_cordinates,alpha,sigma,0);
        %generate BD-users channel gains before SIC order
        G_BD_users_unordered=channelGain_BD(BD_cordinates,users_cordinates,alpha,0);
        
        %%%%% ordering channel gains for the perfect CSI %%%%%
        [G_BS_users,I1]=sort(G_BS_users_unordered,'descend');
        G_BD_users=G_BD_users_unordered(I1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE OF IMPERFECT CSI %%%%%%%%%%%%%%%
        
        %%%%%% channel gains AFFECTED by the error1 %%%%%%%%%
        G_BS_BD_error1=channelGain_BS(0,BD_cordinates,alpha,sigma,error1(2));
        [G_BS_users_error1,B1]=channelGain_BS(1,users_cordinates,alpha,sigma,error1(1));
        G_BD_users_unordered_error1=channelGain_BD(BD_cordinates,users_cordinates,alpha,error1(3));
        G_BD_users_error1=G_BD_users_unordered_error1(B1);
        
        %%%%%% channel gains AFFECTED by the error2 %%%%%%%%%
        G_BS_BD_error2=channelGain_BS(0,BD_cordinates,alpha,sigma,error2(2));
        [G_BS_users_error2,B2]=channelGain_BS(1,users_cordinates,alpha,sigma,error2(1));
        G_BD_users_unordered_error2=channelGain_BD(BD_cordinates,users_cordinates,alpha,error2(3));
        G_BD_users_error2=G_BD_users_unordered_error2(B2);
        
        %%%%%% channel gains AFFECTED by the error3 %%%%%%%%%
        G_BS_BD_error3=channelGain_BS(0,BD_cordinates,alpha,sigma,error3(2));
        [G_BS_users_error3,B3]=channelGain_BS(1,users_cordinates,alpha,sigma,error3(1));
        G_BD_users_unordered_error3=channelGain_BD(BD_cordinates,users_cordinates,alpha,error3(3));
        G_BD_users_error3=G_BD_users_unordered_error3(B3);
        
        %%%%%% channel gains AFFECTED by the error4 %%%%%%%%%
        G_BS_BD_error4=channelGain_BS(0,BD_cordinates,alpha,sigma,error4(2));
        [G_BS_users_error4,B4]=channelGain_BS(1,users_cordinates,alpha,sigma,error4(1));
        G_BD_users_unordered_error4=channelGain_BD(BD_cordinates,users_cordinates,alpha,error4(3));
        G_BD_users_error4=G_BD_users_unordered_error3(B4);
        
        %%%%%% channel gains AFFECTED by the error5 %%%%%%%%%
        G_BS_BD_error5=channelGain_BS(0,BD_cordinates,alpha,sigma,error5(2));
        [G_BS_users_error5,B5]=channelGain_BS(1,users_cordinates,alpha,sigma,error5(1));
        G_BD_users_unordered_error5=channelGain_BD(BD_cordinates,users_cordinates,alpha,error5(3));
        G_BD_users_error5=G_BD_users_unordered_error3(B5);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE OF PERFECT CSI %%%%%%%%%%%%%%%
        
        %calling the function that computes R according to equation (6) in the paper for perfect channel gains
        R = rho_plus(G_BS_users,G_BS_BD,G_BD_users);
        R_unordered = rho_plus(G_BS_users_unordered,G_BS_BD,G_BD_users_unordered);
        
        %compute the optimal rho for ordered channels
        if (isempty(R))
            rho_NOMA=1;
        else
            rho_NOMA=min(1,min(R));
        end
        %compute the optimal rho for unordered channels
        if (isempty(R_unordered))
            rho_NOMA_unordered=1;
        else
            rho_NOMA_unordered=min(1,min(R_unordered));
        end
        
        %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
        
        %conventional NOMA ordered
        G_BS_users=G_BS_users.^2; %channel gain : |h|^2
        Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K(i)));
        end
        
        %conventional NOMA unordered
        G_BS_users_unordered=G_BS_users_unordered.^2; %channel gain : |h|^2
        Pmin_NOMA_conv_unordered=(A(end)-1)/G_BS_users_unordered(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv_unordered=Pmin_NOMA_conv_unordered+(A(j)-1)/G_BS_users_unordered(j)*prod(A((j+1):K(i)));
        end
        
        % NOMA with backscatter ordered
        G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2; %channel gain : |h|^2
        Pmin_NOMA_BD=(A(end)-1)/G_NOMA_BD(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD=Pmin_NOMA_BD+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:K(i)));
        end
        
        % NOMA with backscatter unordered
        G_NOMA_BD_unordered=(sqrt(G_BS_users_unordered)+sqrt(rho_NOMA_unordered)*G_BS_BD*G_BD_users_unordered).^2; %channel gain : |h|^2
        Pmin_NOMA_BD_unordered=(A(end)-1)/G_NOMA_BD_unordered(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD_unordered=Pmin_NOMA_BD_unordered+(A(j)-1)/G_NOMA_BD_unordered(j)*prod(A(j+1:K(i)));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE OF IMPERFECT CSI %%%%%%%%%%%%%%%
        
        %%%%%%% error 1 %%%%%%%%
        
        R_error1= rho_plus(G_BS_users_error1,G_BS_BD_error1,G_BD_users_error1);
        %compute the optimal rho AFTER error1
        if (isempty(R_error1))
            rho_NOMA_error1=1;
        else
            rho_NOMA_error1=min(1,min(R_error1));
        end
        
        %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
        
        %conventional NOMA with error1
        G_BS_users_error1=G_BS_users_error1.^2;%channel gain : |h|^2
        Pmin_NOMA_conv_error1=(A(end)-1)/G_BS_users_error1(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv_error1=Pmin_NOMA_conv_error1+(A(j)-1)/G_BS_users_error1(j)*prod(A((j+1):K(i)));
        end
        
        %NOMA with backscatter with error1
        G_NOMA_BD_error1=(sqrt(G_BS_users_error1)+sqrt(rho_NOMA_error1)*G_BS_BD_error1*G_BD_users_error1).^2;%channel gain : |h|^2
        Pmin_NOMA_BD_error1=(A(end)-1)/G_NOMA_BD_error1(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD_error1=Pmin_NOMA_BD_error1+(A(j)-1)/G_NOMA_BD_error1(j)*prod(A(j+1:K(i)));
        end
        
        
        %%%%%%%%%%%%%%  error 2  %%%%%%%%%%%%%
        
        R_error2= rho_plus(G_BS_users_error2,G_BS_BD_error2,G_BD_users_error2);
        %compute the optimal rho AFTER error2
        if (isempty(R_error2))
            rho_NOMA_error2=1;
        else
            rho_NOMA_error2=min(1,min(R_error2));
        end
        
        %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
        
        %conventional NOMA with error2
        G_BS_users_error2=G_BS_users_error2.^2;%channel gain : |h|^2
        Pmin_NOMA_conv_error2=(A(end)-1)/G_BS_users_error2(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv_error2=Pmin_NOMA_conv_error2+(A(j)-1)/G_BS_users_error2(j)*prod(A((j+1):K(i)));
        end
        
        %NOMA with backscatter with error2
        G_NOMA_BD_error2=(sqrt(G_BS_users_error2)+sqrt(rho_NOMA_error2)*G_BS_BD_error2*G_BD_users_error2).^2;%channel gain : |h|^2
        Pmin_NOMA_BD_error2=(A(end)-1)/G_NOMA_BD_error2(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD_error2=Pmin_NOMA_BD_error2+(A(j)-1)/G_NOMA_BD_error2(j)*prod(A(j+1:K(i)));
        end
        
        %%%%%%%%%%%%%%  error 3  %%%%%%%%%%%%%
        
        R_error3= rho_plus(G_BS_users_error3,G_BS_BD_error3,G_BD_users_error2);
        %compute the optimal rho AFTER error
        if (isempty(R_error3))
            rho_NOMA_error3=1;
        else
            rho_NOMA_error3=min(1,min(R_error3));
        end
        
        %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
        
        %conventional NOMA with error3
        G_BS_users_error3=G_BS_users_error3.^2;%channel gain : |h|^2
        Pmin_NOMA_conv_error3=(A(end)-1)/G_BS_users_error3(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv_error3=Pmin_NOMA_conv_error3+(A(j)-1)/G_BS_users_error3(j)*prod(A((j+1):K(i)));
        end
        
        %NOMA with backscatter with error 3
        G_NOMA_BD_error3=(sqrt(G_BS_users_error3)+sqrt(rho_NOMA_error3)*G_BS_BD_error3*G_BD_users_error3).^2;
        Pmin_NOMA_BD_error3=(A(end)-1)/G_NOMA_BD_error3(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD_error3=Pmin_NOMA_BD_error3+(A(j)-1)/G_NOMA_BD_error3(j)*prod(A(j+1:K(i)));
        end
        
        %%%%%%%%%%%%%%  error 4  %%%%%%%%%%%%%
        
        R_error4= rho_plus(G_BS_users_error4,G_BS_BD_error4,G_BD_users_error4);
        %compute the optimal rho AFTER error
        if (isempty(R_error4))
            rho_NOMA_error4=1;
        else
            rho_NOMA_error4=min(1,min(R_error4));
        end
        
        %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
        
        %conventional NOMA with error4
        G_BS_users_error4=G_BS_users_error4.^2;%channel gain : |h|^2
        Pmin_NOMA_conv_error4=(A(end)-1)/G_BS_users_error4(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv_error4=Pmin_NOMA_conv_error4+(A(j)-1)/G_BS_users_error4(j)*prod(A((j+1):K(i)));
        end
        
        %NOMA with backscatter with error 4
        G_NOMA_BD_error4=(sqrt(G_BS_users_error4)+sqrt(rho_NOMA_error4)*G_BS_BD_error4*G_BD_users_error4).^2;
        Pmin_NOMA_BD_error4=(A(end)-1)/G_NOMA_BD_error4(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD_error4=Pmin_NOMA_BD_error4+(A(j)-1)/G_NOMA_BD_error4(j)*prod(A(j+1:K(i)));
        end
        
        
        %%%%%%%%%%%%%%  error 5  %%%%%%%%%%%%%
        
        R_error5= rho_plus(G_BS_users_error5,G_BS_BD_error5,G_BD_users_error5);
        %compute the optimal rho AFTER error
        if (isempty(R_error5))
            rho_NOMA_error5=1;
        else
            rho_NOMA_error5=min(1,min(R_error5));
        end
        
        %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
        
        %conventional NOMA with error5
        G_BS_users_error5=G_BS_users_error5.^2;%channel gain : |h|^2
        Pmin_NOMA_conv_error5=(A(end)-1)/G_BS_users_error5(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_conv_error5=Pmin_NOMA_conv_error5+(A(j)-1)/G_BS_users_error5(j)*prod(A((j+1):K(i)));
        end
        
        %NOMA with backscatter with error 5
        G_NOMA_BD_error5=(sqrt(G_BS_users_error5)+sqrt(rho_NOMA_error5)*G_BS_BD_error5*G_BD_users_error5).^2;
        Pmin_NOMA_BD_error5=(A(end)-1)/G_NOMA_BD_error5(end);
        for j=1:(K(i)-1)
            Pmin_NOMA_BD_error5=Pmin_NOMA_BD_error5+(A(j)-1)/G_NOMA_BD_error5(j)*prod(A(j+1:K(i)));
        end
        
        
        %---- checking the feasability condition ----%
        while (Pmin_NOMA_conv>Pmax || Pmin_NOMA_conv_unordered >Pmax || Pmin_NOMA_BD>Pmax ||Pmin_NOMA_BD_unordered>Pmax || Pmin_NOMA_conv_error1>Pmax || Pmin_NOMA_conv_error2>Pmax || Pmin_NOMA_conv_error3>Pmax || Pmin_NOMA_conv_error4>Pmax || Pmin_NOMA_conv_error5>Pmax || Pmin_NOMA_BD_error1>Pmax || Pmin_NOMA_BD_error2>Pmax || Pmin_NOMA_BD_error3>Pmax|| Pmin_NOMA_BD_error4>Pmax || Pmin_NOMA_BD_error5>Pmax)
            
            %---- generating channel gains ----%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE OF PERFECT CSI %%%%%%%%%%%%%%%

            %generate x and y coordinates for users
            users_cordinates = coordinates(K(i),radius_BS_users,min_dis_BS_users)';
            %generate x and y coordinates for BD
            BD_cordinates = coordinates(1,radius_BS_BD,min_dis_BS_BD);
            
            %generate BS-BD channel gain
            G_BS_BD=channelGain_BS(0,BD_cordinates,alpha,sigma,0);
            %generate BS-users channel gains (in descending order -> SIC)
            [G_BS_users_unordered,I]=channelGain_BS(0,users_cordinates,alpha,sigma,0);
            %generate BD-users channel gains before SIC order
            G_BD_users_unordered=channelGain_BD(BD_cordinates,users_cordinates,alpha,0);
            
            %%%%% ordering channel gains for the perfect CSI %%%%%
            [G_BS_users,I1]=sort(G_BS_users_unordered,'descend');
            G_BD_users=G_BD_users_unordered(I1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE OF IMPERFECT CSI %%%%%%%%%%%%%%%
            
            %%%%%% channel gains AFFECTED by the error1 %%%%%%%%%
            G_BS_BD_error1=channelGain_BS(0,BD_cordinates,alpha,sigma,error1(2));
            [G_BS_users_error1,B1]=channelGain_BS(1,users_cordinates,alpha,sigma,error1(1));
            G_BD_users_unordered_error1=channelGain_BD(BD_cordinates,users_cordinates,alpha,error1(3));
            G_BD_users_error1=G_BD_users_unordered_error1(B1);
            
            %%%%%% channel gains AFFECTED by the error2 %%%%%%%%%
            G_BS_BD_error2=channelGain_BS(0,BD_cordinates,alpha,sigma,error2(2));
            [G_BS_users_error2,B2]=channelGain_BS(1,users_cordinates,alpha,sigma,error2(1));
            G_BD_users_unordered_error2=channelGain_BD(BD_cordinates,users_cordinates,alpha,error2(3));
            G_BD_users_error2=G_BD_users_unordered_error2(B2);
            
            %%%%%% channel gains AFFECTED by the error3 %%%%%%%%%
            G_BS_BD_error3=channelGain_BS(0,BD_cordinates,alpha,sigma,error3(2));
            [G_BS_users_error3,B3]=channelGain_BS(1,users_cordinates,alpha,sigma,error3(1));
            G_BD_users_unordered_error3=channelGain_BD(BD_cordinates,users_cordinates,alpha,error3(3));
            G_BD_users_error3=G_BD_users_unordered_error3(B3);
            
            %%%%%% channel gains AFFECTED by the error4 %%%%%%%%%
            G_BS_BD_error4=channelGain_BS(0,BD_cordinates,alpha,sigma,error4(2));
            [G_BS_users_error4,B4]=channelGain_BS(1,users_cordinates,alpha,sigma,error4(1));
            G_BD_users_unordered_error4=channelGain_BD(BD_cordinates,users_cordinates,alpha,error4(3));
            G_BD_users_error4=G_BD_users_unordered_error3(B4);
            
            %%%%%% channel gains AFFECTED by the error5 %%%%%%%%%
            G_BS_BD_error5=channelGain_BS(0,BD_cordinates,alpha,sigma,error5(2));
            [G_BS_users_error5,B5]=channelGain_BS(1,users_cordinates,alpha,sigma,error5(1));
            G_BD_users_unordered_error5=channelGain_BD(BD_cordinates,users_cordinates,alpha,error5(3));
            G_BD_users_error5=G_BD_users_unordered_error3(B5);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE OF PERFECT CSI %%%%%%%%%%%%%%%
            %calling the function that computes R according to equation (6) in the paper for perfect channel gains
            R = rho_plus(G_BS_users,G_BS_BD,G_BD_users);
            R_unordered = rho_plus(G_BS_users_unordered,G_BS_BD,G_BD_users_unordered);
            
            %compute the optimal rho for ordered channels
            if (isempty(R))
                rho_NOMA=1;
            else
                rho_NOMA=min(1,min(R));
            end
            %compute the optimal rho for unordered channels
            if (isempty(R_unordered))
                rho_NOMA_unordered=1;
            else
                rho_NOMA_unordered=min(1,min(R_unordered));
            end
            
            %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
            
            %conventional NOMA ordered
            G_BS_users=G_BS_users.^2; %channel gain : |h|^2
            Pmin_NOMA_conv=(A(end)-1)/G_BS_users(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv=Pmin_NOMA_conv+(A(j)-1)/G_BS_users(j)*prod(A((j+1):K(i)));
            end
            
            %conventional NOMA unordered
            G_BS_users_unordered=G_BS_users_unordered.^2; %channel gain : |h|^2
            Pmin_NOMA_conv_unordered=(A(end)-1)/G_BS_users_unordered(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv_unordered=Pmin_NOMA_conv_unordered+(A(j)-1)/G_BS_users_unordered(j)*prod(A((j+1):K(i)));
            end
            
            % NOMA with backscatter ordered
            G_NOMA_BD=(sqrt(G_BS_users)+sqrt(rho_NOMA)*G_BS_BD*G_BD_users).^2; %channel gain : |h|^2
            Pmin_NOMA_BD=(A(end)-1)/G_NOMA_BD(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD=Pmin_NOMA_BD+(A(j)-1)/G_NOMA_BD(j)*prod(A(j+1:K(i)));
            end
            
            % NOMA with backscatter unordered
            G_NOMA_BD_unordered=(sqrt(G_BS_users_unordered)+sqrt(rho_NOMA_unordered)*G_BS_BD*G_BD_users_unordered).^2; %channel gain : |h|^2
            Pmin_NOMA_BD_unordered=(A(end)-1)/G_NOMA_BD_unordered(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD_unordered=Pmin_NOMA_BD_unordered+(A(j)-1)/G_NOMA_BD_unordered(j)*prod(A(j+1:K(i)));
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE OF IMPERFECT CSI %%%%%%%%%%%%%%%
            
            %%%%%%% error 1 %%%%%%%%
            
            R_error1= rho_plus(G_BS_users_error1,G_BS_BD_error1,G_BD_users_error1);
            %compute the optimal rho AFTER error1
            if (isempty(R_error1))
                rho_NOMA_error1=1;
            else
                rho_NOMA_error1=min(1,min(R_error1));
            end
            
            %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
            
            %conventional NOMA with error1
            G_BS_users_error1=G_BS_users_error1.^2;%channel gain : |h|^2
            Pmin_NOMA_conv_error1=(A(end)-1)/G_BS_users_error1(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv_error1=Pmin_NOMA_conv_error1+(A(j)-1)/G_BS_users_error1(j)*prod(A((j+1):K(i)));
            end
            
            %NOMA with backscatter with error1
            G_NOMA_BD_error1=(sqrt(G_BS_users_error1)+sqrt(rho_NOMA_error1)*G_BS_BD_error1*G_BD_users_error1).^2;%channel gain : |h|^2
            Pmin_NOMA_BD_error1=(A(end)-1)/G_NOMA_BD_error1(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD_error1=Pmin_NOMA_BD_error1+(A(j)-1)/G_NOMA_BD_error1(j)*prod(A(j+1:K(i)));
            end
            
            
            %%%%%%%%%%%%%%  error 2  %%%%%%%%%%%%%
            
            R_error2= rho_plus(G_BS_users_error2,G_BS_BD_error2,G_BD_users_error2);
            %compute the optimal rho AFTER error2
            if (isempty(R_error2))
                rho_NOMA_error2=1;
            else
                rho_NOMA_error2=min(1,min(R_error2));
            end
            
            %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
            
            %conventional NOMA with error2
            G_BS_users_error2=G_BS_users_error2.^2;%channel gain : |h|^2
            Pmin_NOMA_conv_error2=(A(end)-1)/G_BS_users_error2(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv_error2=Pmin_NOMA_conv_error2+(A(j)-1)/G_BS_users_error2(j)*prod(A((j+1):K(i)));
            end
            
            %NOMA with backscatter with error2
            G_NOMA_BD_error2=(sqrt(G_BS_users_error2)+sqrt(rho_NOMA_error2)*G_BS_BD_error2*G_BD_users_error2).^2;%channel gain : |h|^2
            Pmin_NOMA_BD_error2=(A(end)-1)/G_NOMA_BD_error2(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD_error2=Pmin_NOMA_BD_error2+(A(j)-1)/G_NOMA_BD_error2(j)*prod(A(j+1:K(i)));
            end
            
            %%%%%%%%%%%%%%  error 3  %%%%%%%%%%%%%
            
            R_error3= rho_plus(G_BS_users_error3,G_BS_BD_error3,G_BD_users_error2);
            %compute the optimal rho AFTER error
            if (isempty(R_error3))
                rho_NOMA_error3=1;
            else
                rho_NOMA_error3=min(1,min(R_error3));
            end
            
            %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
            
            %conventional NOMA with error3
            G_BS_users_error3=G_BS_users_error3.^2;%channel gain : |h|^2
            Pmin_NOMA_conv_error3=(A(end)-1)/G_BS_users_error3(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv_error3=Pmin_NOMA_conv_error3+(A(j)-1)/G_BS_users_error3(j)*prod(A((j+1):K(i)));
            end
            
            %NOMA with backscatter with error 3
            G_NOMA_BD_error3=(sqrt(G_BS_users_error3)+sqrt(rho_NOMA_error3)*G_BS_BD_error3*G_BD_users_error3).^2;
            Pmin_NOMA_BD_error3=(A(end)-1)/G_NOMA_BD_error3(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD_error3=Pmin_NOMA_BD_error3+(A(j)-1)/G_NOMA_BD_error3(j)*prod(A(j+1:K(i)));
            end
            
            %%%%%%%%%%%%%%  error 4  %%%%%%%%%%%%%
            
            R_error4= rho_plus(G_BS_users_error4,G_BS_BD_error4,G_BD_users_error4);
            %compute the optimal rho AFTER error
            if (isempty(R_error4))
                rho_NOMA_error4=1;
            else
                rho_NOMA_error4=min(1,min(R_error4));
            end
            
            %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
            
            %conventional NOMA with error4
            G_BS_users_error4=G_BS_users_error4.^2;%channel gain : |h|^2
            Pmin_NOMA_conv_error4=(A(end)-1)/G_BS_users_error4(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv_error4=Pmin_NOMA_conv_error4+(A(j)-1)/G_BS_users_error4(j)*prod(A((j+1):K(i)));
            end
            
            %NOMA with backscatter with error 4
            G_NOMA_BD_error4=(sqrt(G_BS_users_error4)+sqrt(rho_NOMA_error4)*G_BS_BD_error4*G_BD_users_error4).^2;
            Pmin_NOMA_BD_error4=(A(end)-1)/G_NOMA_BD_error4(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD_error4=Pmin_NOMA_BD_error4+(A(j)-1)/G_NOMA_BD_error4(j)*prod(A(j+1:K(i)));
            end
            
            
            %%%%%%%%%%%%%%  error 5  %%%%%%%%%%%%%
            
            R_error5= rho_plus(G_BS_users_error5,G_BS_BD_error5,G_BD_users_error5);
            %compute the optimal rho AFTER error
            if (isempty(R_error5))
                rho_NOMA_error5=1;
            else
                rho_NOMA_error5=min(1,min(R_error5));
            end
            
            %---compute the minimum power budget (Pmin) required for meeting QoS constraints---%
            
            %conventional NOMA with error5
            G_BS_users_error5=G_BS_users_error5.^2;%channel gain : |h|^2
            Pmin_NOMA_conv_error5=(A(end)-1)/G_BS_users_error5(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_conv_error5=Pmin_NOMA_conv_error5+(A(j)-1)/G_BS_users_error5(j)*prod(A((j+1):K(i)));
            end
            
            %NOMA with backscatter with error 5
            G_NOMA_BD_error5=(sqrt(G_BS_users_error5)+sqrt(rho_NOMA_error5)*G_BS_BD_error5*G_BD_users_error5).^2;
            Pmin_NOMA_BD_error5=(A(end)-1)/G_NOMA_BD_error5(end);
            for j=1:(K(i)-1)
                Pmin_NOMA_BD_error5=Pmin_NOMA_BD_error5+(A(j)-1)/G_NOMA_BD_error5(j)*prod(A(j+1:K(i)));
            end
            
        end
        
        %---- compute the optimal energy efficiency for PERFTECT CSI ----%
        [EE_NOMA_BD(i),outage_perfectCSI(i)] = Dinkelbach_NOMA_outage(q,G_BS_users,G_NOMA_BD,A,Pmax,Pmin_NOMA_conv,Pc);
        
        %---- compute the optimal energy efficiency and corresponding outage for IMPERFTECT CSI ----%
        [EE_NOMA_BD_error1(i),outage1(i)] = Dinkelbach_NOMA_ICSI(q,G_BS_users,G_NOMA_BD,G_BS_users_error1,G_NOMA_BD_error1,A,Pmax,Pmin_NOMA_conv_error1,Pc);
        [EE_NOMA_BD_error2(i),outage2(i)] = Dinkelbach_NOMA_ICSI(q,G_BS_users,G_NOMA_BD,G_BS_users_error2,G_NOMA_BD_error2,A,Pmax,Pmin_NOMA_conv_error2,Pc);
        [EE_NOMA_BD_error3(i),outage3(i)] = Dinkelbach_NOMA_ICSI(q,G_BS_users,G_NOMA_BD,G_BS_users_error3,G_NOMA_BD_error3,A,Pmax,Pmin_NOMA_conv_error3,Pc);
        [EE_NOMA_BD_error4(i),outage4(i)] = Dinkelbach_NOMA_ICSI(q,G_BS_users,G_NOMA_BD,G_BS_users_error4,G_NOMA_BD_error4,A,Pmax,Pmin_NOMA_conv_error4,Pc);
        [EE_NOMA_BD_error5(i),outage5(i)] = Dinkelbach_NOMA_ICSI(q,G_BS_users,G_NOMA_BD,G_BS_users_error5,G_NOMA_BD_error5,A,Pmax,Pmin_NOMA_conv_error5,Pc);
        
    end
    %-------- stock results for each channel realization --------%
    EE_NOMA_BD_n(n,:)=EE_NOMA_BD;
    EE_NOMA_BD_error1_n(n,:)=EE_NOMA_BD_error1;
    EE_NOMA_BD_error2_n(n,:)=EE_NOMA_BD_error2;
    EE_NOMA_BD_error3_n(n,:)=EE_NOMA_BD_error3;
    EE_NOMA_BD_error4_n(n,:)=EE_NOMA_BD_error4;
    EE_NOMA_BD_error5_n(n,:)=EE_NOMA_BD_error5;
    
    outage_perfectCSI_n(n,:)=outage_perfectCSI;
    outage1_n(n,:)=outage1;
    outage2_n(n,:)=outage2;
    outage3_n(n,:)=outage3;
    outage4_n(n,:)=outage4;
    outage5_n(n,:)=outage5;
    
end

% summation of outage
sum_outage1=sum(outage1_n,1);
sum_outage2=sum(outage2_n,1);
sum_outage3=sum(outage3_n,1);
sum_outage4=sum(outage4_n,1);
sum_outage5=sum(outage5_n,1);

%%%% to avoid dividing by 0
for i=1:length(K)
    if (N==sum_outage1(i))
    count1(i)=1;
    else
        count1(i)=N-sum_outage1(i);
    end
end

for i=1:length(K)
    if (N==sum_outage2(i))
    count2(i)=1;
    else
        count2(i)=N-sum_outage2(i);
    end
end

for i=1:length(K)
    if (N==sum_outage3(i))
    count3(i)=1;
    else
        count3(i)=N-sum_outage3(i);
    end
end

for i=1:length(K)
    if (N==sum_outage4(i))
    count4(i)=1;
    else
        count4(i)=N-sum_outage4(i);
    end
end

for i=1:length(K)
    if (N==sum_outage5(i))
    count5(i)=1;
    else
        count5(i)=N-sum_outage5(i);
    end
end


%------------- averaging over channel realizations -----------%
EE_NOMA_BD_mean=sum(EE_NOMA_BD_n)./(N-sum(outage_perfectCSI_n,1));
EE_NOMA_BD_error1_mean=sum(EE_NOMA_BD_error1_n,1)./count1;
EE_NOMA_BD_error2_mean=sum(EE_NOMA_BD_error2_n,1)./count2;
EE_NOMA_BD_error3_mean=sum(EE_NOMA_BD_error3_n,1)./count3;
EE_NOMA_BD_error4_mean=sum(EE_NOMA_BD_error4_n,1)./N;
EE_NOMA_BD_error5_mean=sum(EE_NOMA_BD_error5_n,1)./N;

%%%% computing the empirical outage
outage_perfectCSI_mean=sum(outage_perfectCSI_n)./N;
outage1_mean=100*sum(outage1_n)./N;
outage2_mean=100*sum(outage2_n)./N;
outage3_mean=100*sum(outage3_n)./N;
outage4_mean=100*sum(outage4_n)./N;
outage5_mean=100*sum(outage5_n)./N;


%------------- plot figures -----------%
figure(1)
plot(K,EE_NOMA_BD_mean,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_NOMA_BD_error1_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_NOMA_BD_error2_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_NOMA_BD_error3_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_NOMA_BD_error4_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,EE_NOMA_BD_error5_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
ylabel('\xi_{EE} (bits/J)');
xlabel('Number of receivers');
grid on;
legend('perfect CSI', 'imperfect CSI \sigma_h^2/\sigma_e^2=20dB','imperfect CSI \sigma_h^2/\sigma_e^2=10dB','imperfect CSI \sigma_h^2/\sigma_e^2=0dB','imperfect CSI \sigma_h^2/\sigma_e^2=-10dB','imperfect CSI \sigma_h^2/\sigma_e^2=-20dB','Location','NorthEast');
%a=axes('position',get(gca,'position'),'visible','off');
%legend(a, [d4 d5 d6],'imperfect CSI \sigma_h^2/\sigma_e^2=0dB','imperfect CSI \sigma_h^2/\sigma_e^2=-10dB','imperfect CSI \sigma_h^2/\sigma_e^2=-20dB','Location','East');


figure(2)
plot(K,outage_perfectCSI_mean,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,outage1_mean,'-o','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,outage2_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,outage3_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,outage4_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
hold on;
plot(K,outage5_mean,'->','MarkerSize',7,'LineWidth',3,'MarkerIndices', 1:length(K));
ylabel('Outage (%)');
xlabel('Number of receivers');
grid on;
legend('perfect CSI','imperfect CSI \sigma_h^2/\sigma_e^2=20dB', 'imperfect CSI \sigma_h^2/\sigma_e^2=10dB','imperfect CSI \sigma_h^2/\sigma_e^2=0dB','imperfect CSI \sigma_h^2/\sigma_e^2=-10dB','imperfect CSI \sigma_h^2/\sigma_e^2=-20dB','Location','best');
grid on;
