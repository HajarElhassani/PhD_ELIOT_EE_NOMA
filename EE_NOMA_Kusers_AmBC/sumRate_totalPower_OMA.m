%% This function computes the sum rate, total power consumption and the energy efficiency as a tradeoff between the two for a certain alpha
% input:          H -> channel gain of the direct link
%                 G -> channel gain of the backscattered+direct link
%                 q -> probability of backscattering (B=1)
%                 A -V A=2^(2*Rmin) the SNR threshold
%                 Pmax -> maximum power budget
%                 Pc -> circuit power (Pc)
%                 alpha -> the tradeoff parameter
% ouput:          EE -> energy efficiency as a tradeoff
%                 sumRate -> the sum rate
%                 totalPower -> the total power consumption 

function [sumRate,totalPower,EE] = sumRate_totalPower_OMA(q,alpha,H,G,A,Pmax,Pc)

% find the optimal power allocation for certain alpha
p = compute_EE_OMA(q,alpha,Pmax,A,H,G,Pc);
%compute the rate
r=(q/(2*length(G)))*log2(1+G.*p')+((1-q)/(2*length(H)))*log2(1+H.*p');
%compute the sum rate
sumRate=sum(r);
% compute the total power consumption
totalPower=10*log10((sum(p)/length(H)+Pc))+30;
%compute the energy efficiency as a tradeoff
EE=sumRate-alpha*totalPower;

end
