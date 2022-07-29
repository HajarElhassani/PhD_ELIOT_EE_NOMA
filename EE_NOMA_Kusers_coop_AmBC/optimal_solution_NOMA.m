%% This function compute the GEE according to Dinkelbach's algorithm

%% input: channel gain (G), A=2^(2*Rmin), maximum power budget (Pmax), circuit power (Pc)
%% ouput: global energy efficiency (alpha)


function alpha = optimal_solution_NOMA(G,A,Pmax,Pmin,Pc)
%% initialization parameters
epsilon=1e-7;
alpha=1e-5;
R_P=1;
%% lower and upper bound of feasible p1
lb=(A(1)-1)/G(1); %lower bound
ub=1/prod(A(2:length(G)))*(Pmax-Pmin+lb*prod(A(2:length(G)))); % upper bound

%% Dinkelbach's algorithm 
while(R_P>epsilon)
    c=1/(2*log(2)*alpha*prod(A(2:length(G))))-1/G(1);% critical point
    %% compute the optimal power allocation
    theta(1)=min(ub,max(c,lb));
    theta(2)=(A(2)-1)/G(2)+theta(1)*A(2);
    
    if length(G)>2
        for k=3:length(G)
            B=0;
            for i=2:k-1
                B=B+(A(i)-1)/G(i)*prod(A(i+1:k));
            end
            theta(k)=(A(k)-1)/G(k)+theta(1)*prod(A(2:k))+B;
        end
    end
    %% compute the rates
    r=zeros(length(G),1);
    r(1)=1/2*log2(1+G(1)*theta(1));
    for k=2:length(G)
        r(k)=1/2*log2((1+G(k)*theta(k))/(1+G(k)*theta(k-1)));
    end
    %% compute the function sum rate vs total power tradeoff
    R_P=sum(r)-alpha*(theta(end)+Pc);  
    %% update alpha
    alpha=sum(r)/(theta(end)+Pc);

end

end