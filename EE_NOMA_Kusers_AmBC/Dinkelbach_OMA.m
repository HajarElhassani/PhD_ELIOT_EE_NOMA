%------ Dinkelbach's algorithm -------% 

% input:          H -> channel gain of the direct link
%                 G -> channel gain of the backscattered+direct link
%                 A -V A=2^(2*Rmin) the SNR threshold
%                 Pmax -> maximum power budget 
%                 Pc -> circuit power (Pc)
% ouput:         alpha -> global energy efficiency (alpha)


function [alpha] = Dinkelbach_OMA(q,H,G,A,Pmax,Pc)

%stop criterion for Dinkelbach's algorithm
epsilon=1e-7;
%initilization of alpha
alpha=1e-5;
%initialization of Energy Efficiency as tradeoff between sum rate and power consumption
R_P=1;

while(R_P>epsilon)
    % compute the optimal power
    p=compute_EE_OMA(q,alpha,Pmax,A,H,G,Pc);    
    % compute the rates
    r=(q/(2*length(G)))*log2(1+G.*p')+((1-q)/(2*length(H)))*log2(1+H.*p');
    % compute the sum rate vs total power consumption tradeoff
    R_P=sum(r)-alpha*(sum(p)/length(H)+Pc);
    % compute alpha (energy efficiency)
    alpha=sum(r)/(sum(p)/length(H)+Pc);
end

end