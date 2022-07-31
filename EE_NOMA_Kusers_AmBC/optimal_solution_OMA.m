%% This function computes the energy efficiency in the case of q=0 (without backscattered which has simpler solution) using dinkelbach's algorithm

% input            H -> the channel gain of the direct link 
%                  A -> SNR threshold =2^(2*Rmin)
%                  Pmax -> maximum power budget
%                  Pc -> circuit power
% ouput            alpha -> energy efficiency
%                  Rsum -> the sum rate
%                  Ptot -> the total power consumption


function [alpha,Rsum,Ptot] = optimal_solution_OMA(H,A,Pmax,Pc)

% initialization
epsilon=1e-7;
alpha=1e-5;
R_P=1;

% lower and upper bound of the power allocation
P1=(A.^length(H)-1)./H; % lower bound
P3=Pmax*ones(length(H),1); %upper bound

while(R_P>epsilon)

    P2=1/(2*log(2)*alpha)-1./H;% critical point

    for k=1:length(H)
        p(k)=min(P3(k),max(P1(k),P2(k)));
    end
    %compute the sum rate vs power consumption tradeoff
    R_P=(1/(2*length(H)))*sum(log2(1+H.*p'))-alpha*(sum(p)/length(H)+Pc);
    %update alpha(energy efficiency)
    alpha=(1/(2*length(H)))*sum(log2(1+H.*p'))/(sum(p)/length(H)+Pc);
    
    %the sum rate
    Rsum=(1/(2*length(H)))*sum(log2(1+H.*p'));
    % the total power consumption
    Ptot=(sum(p)/length(H))+Pc;
end
end