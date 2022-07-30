%% this function computes NOMA energy efficiency ratio using Dinkelbach's algorithm
%input: G: channel gain, A=2^(2Rmin), Pmax=power budget, Pc: circuit power
%output: alpha: energy efficiency
function alpha = optimal_solution_NOMA(G,A,Pmax,Pmin,Pc)

%initialization
epsilon=1e-7;
alpha=1e-5;
R_P=1;

%lower bound and upper bounds on p1 
lb=(A(1)-1)/G(1);
ub=1/prod(A(2:length(G)))*(Pmax-Pmin+lb*prod(A(2:length(G))));

while(R_P>epsilon)
    
    %compute the power allocation through \theta(k)
    criticalPoint=1/(2*log(2)*alpha*prod(A(2:length(G))))-1/G(1);
    theta(1)=min(ub,max(criticalPoint,lb));
    theta(2)=(A(2)-1)/G(2)+theta(1)*A(2);
    if length(G)>2
        for k=3:length(G)
            B=0;
            for i=2:k-1
                B=B+(A(i)-1)/G(i)*prod(A(i+1:k));
            end
            theta(k)=(A(k)-1)/G(k)+B+theta(1)*prod(A(2:k));
        end
    end
    %initialize the rates
    r=zeros(length(G),1);
    %compute the rates
    r(1)=1/2*log2(1+G(1)*theta(1));
    for k=2:length(G)
        r(k)=1/2*log2((1+G(k)*theta(k))/(1+G(k)*theta(k-1)));
    end
    
    %compute F(alpha)=sumRate-alpha*powerConsumption
    R_P=sum(r)-alpha*(theta(end)+Pc);
    alpha=sum(r)/(theta(end)+Pc);
    
    %until convergence
end

end
