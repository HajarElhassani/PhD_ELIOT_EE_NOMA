%------ Dinkelbach's algorithm -------% 

% input:          H -> channel gain of the direct link
%                 G -> channel gain of the backscattered+direct link
%                 A -V A=2^(2*Rmin) the SNR threshold
%                 Pmax -> maximum power budget 
%                 Pc -> circuit power (Pc)
%                 q -> the probability of backscattering (B=1)
% ouput:         alpha -> global energy efficiency (alpha)

function [alpha] = Dinkelbach_NOMA(q,H,G,A,Pmax,Pmin,Pc)

%---- initialization ----%
epsilon=1e-7;
alpha=1e-5;
R_P=1;
r=zeros(1,length(H));

% the upper and lower bound on the feasible p1
lb=(A(1)-1)/H(1); %lower bound
ub=1/prod(A(2:length(H)))*(Pmax-Pmin+lb*prod(A(2:length(H)))); %upper bound

%----- Dinkelbach's algorithm ----%
while(R_P>epsilon)
    %compute the optimal power p1 and the corresponding energy efficiency
    p1=compute_EE_NOMA(q,alpha,A,H,G,Pc,ub,lb);
    %compute the optimal power allocation
    theta(1)=p1;
    theta(2)=(A(2)-1)/H(2)+theta(1)*A(2);
    
    if length(H)>2
        for k=3:length(H)
            B=0;
            for i=2:k-1
                B=B+(A(i)-1)/H(i)*prod(A(i+1:k));
            end
            theta(k)=(A(k)-1)/H(k)+B+theta(1)*prod(A(2:k));
        end
    end
   
    %compute the sum rate
    r(1)=q/2*log2(1+G(1)*theta(1))+(1-q)/2*log2(1+H(1)*theta(1));
    for k=2:length(H)
        r(k)=q/2*log2((1+G(k)*theta(k))/(1+G(k)*theta(k-1)))+(1-q)/2*log2((1+H(k)*theta(k))/(1+H(k)*theta(k-1)));
    end
    %compute the tradeoff sum rate vs. total power consumption
    R_P=sum(r)-alpha*(theta(end)+Pc); 
    %update alpha (energy efficiency)
    alpha=sum(r)/(theta(end)+Pc);

end

end
