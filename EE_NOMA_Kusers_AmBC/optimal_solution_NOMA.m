%% This function computes the energy efficiency in the case of q=0 (without backscattered which has simpler solution) using dinkelbach's algorithm

% input            H -> the channel gain of the direct link 
%                  A -> SNR threshold =2^(2*Rmin)
%                  Pmax -> maximum power budget
%                  Pc -> circuit power
%                  Pmin -> the minimum power requirement
% ouput            alpha -> energy efficiency 

function alpha = optimal_solution_NOMA(H,A,Pmax,Pmin,Pc)
%intilization
epsilon=1e-7;
alpha=1e-5;
R_P=1;

% lower and upper bound on the feasible p1
lb=(A(1)-1)/H(1); % lower bound
ub=1/prod(A(2:length(H)))*(Pmax-Pmin+lb*prod(A(2:length(H)))); %upper bound

while(R_P>epsilon)
    c=1/(2*log(2)*alpha*prod(A(2:length(H))))-1/H(1); %critical point
    
    % compute the optimal power allocation
    theta(1)=min(ub,max(c,lb));
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
    r=zeros(length(H),1);
    r(1)=1/2*log2(1+H(1)*theta(1));
    for k=2:length(H)
        r(k)=1/2*log2((1+H(k)*theta(k))/(1+H(k)*theta(k-1)));
    end
    %compute the tradeoff between sum rate and total power consumption
    R_P=sum(r)-alpha*(theta(end)+Pc);
    %update alpha (energy efficiency)
    alpha=sum(r)/(theta(end)+Pc);

end

end