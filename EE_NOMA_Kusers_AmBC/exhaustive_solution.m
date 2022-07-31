%% This function computes the optimal energy efficiency for different values of q via exhaustive search
% input           q -> the probability of backscattering (B=1)
%                 H -> the channel gain of the direct link
%                 G -> the channel gain of the direct+backscattered link
%                 Pmax -> power budget at BS
%                 Pc -> circuit power
%                 Rmin -> the rate threshold 
%output           EE -> energy efficiency

function [EE] = exhaustive_solution(q,H,G,Pmax,Pc,Rmin)

% the variables range -> the reflection coefficient rho, and the power allocation
% p1,p2
[rho,P1,P2]=ndgrid(0:1/100:1,0:Pmax/1000:Pmax,0:Pmax/1000:Pmax);

% simplifying expressions
A1=(P1.*(H(1)+sqrt(rho)*G(1)).^2);
B1=(P1.*(H(2)+sqrt(rho)*G(2)).^2);

A2=(P2.*(H(2)+sqrt(rho)*G(2)).^2);
B2=(P2.*(H(1)+sqrt(rho)*G(1)).^2);

D1=(A2./(1+B1));
D2=(B2./(1+A1));

F1=(P2*H(2)^2)./(1+P1*H(2)^2);
F2=(P2*H(1)^2)./(1+P1*H(1)^2);

% the constraints
c1 = (P1+P2)<=Pmax;
c2 = (q/2*log2(1+A1)+(1-q)/2*log2(1+P1.*H(1)^2))>=Rmin;
c3 = (q/2*log2(1+D1)+(1-q)/2*log2(1+F1))>=Rmin;
c4 = (q/2*log2(1+D2)+(1-q)/2*log2(1+F2))>=(q/2*log2(1+D1)+(1-q)/2*log2(1+F1));

% filter the range of variables that meet the constraints
pass=(c1) & (c2) & (c3) & (c4);
p1=P1(pass);
p2=P2(pass);
rho2=rho(pass);

%compute the energy efficiency for those filtered variables by grid search
fun=(q/2*(log2(1+p1.*(H(1)+sqrt(rho2)*G(1)).^2)+log2(1+(p2.*(H(2)+sqrt(rho2)*G(2)).^2)./(1+p1.*(H(2)+sqrt(rho2)*G(2)).^2)))+(1-q)/2*(log2(1+p1*H(1)^2)+log2(1+(p2*H(2)^2)./(1+p1*H(2)^2))))./(p1+p2+Pc); 
%take the maximum (optimal) value of energy efficiency
EE=max(fun(:));

end
