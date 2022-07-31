%% This function computes the sum rate, total power consumption and the energy efficiency as a tradeoff between the two for a certain alpha
% input:          H -> channel gain of the direct link
%                 G -> channel gain of the backscattered+direct link
%                 q -> probability of backscattering (B=1)
%                 A -V A=2^(2*Rmin) the SNR threshold
%                 Pmax -> maximum power budget
%                 Pmin -> the minimum power requirement
%                 Pc -> circuit power (Pc)
%                 alpha -> the tradeoff parameter
% ouput:          EE -> energy efficiency as a tradeoff
%                 sumRate -> the sum rate
%                 totalPower -> the total power consumption 


function [sumRate,totalPower,EE] = sumRate_totalPower_NOMA(q,alpha,H,G,A,Pmax,Pmin,Pc)
%---- initialization ----%
r=zeros(1,length(H));
%lower and upper bounds on the power p1
lb=(A(1)-1)/H(1);
ub=1/prod(A(2:length(H)))*(Pmax-Pmin+lb*prod(A(2:length(H))));
%compute the optimal power p1 and the corresponding energy efficiency
p1=compute_EE_NOMA(q,alpha,A,H,G,Pc,ub,lb);

%compute the optimal power allocation
p(1)=p1;
p(2)=(A(2)-1)*(1/H(2)+p(1));

if length(H)>2
    for k=3:length(H)
        B=0;
        for i=2:k-1
            B=B+(A(i)-1)/H(i)*prod(A(i+1:k-1));
        end
        p(k)=(A(k)-1)*(1/H(k)+p(1)*prod(A(2:k-1))+B);
    end
end

%compute the rates
r(1)=q/2*log2(1+G(1)*p(1))+(1-q)/2*log2(1+H(1)*p(1));
for k=2:length(H)
    r(k)=q/2*log2(1+G(k)*p(k)/(1+G(k)*(sum(p(1:k-1)))))+(1-q)/2*log2(1+H(k)*p(k)/(1+H(k)*(sum(p(1:k-1)))));
end

%compute the sum rate
sumRate=sum(r);
%compute the total power consumption in dBm
totalPower=10*log10((sum(p)+Pc))+30;
% the energy efficiency as a tradeoff
EE=sumRate-alpha*totalPower;

end
