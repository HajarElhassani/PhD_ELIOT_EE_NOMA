%% this function computes OMA energy efficiency ratio using Dinkelbach's algorithm
%input: G: channel gain, A=2^(2Rmin), Pmax=power budget, Pc: circuit power
%output: alpha: energy efficiency
function alpha = optimal_solution_OMA(G,A,Pmax,Pc)

%initialization
epsilon=1e-7;
alpha=1e-5;
R_P=1;

%lower bound and upper bound on the powers
ub=(A.^length(G)-1)./G;
lb=Pmax*ones(length(G),1);

%start of Dinkelbach's algorithm
while(R_P>epsilon)
    criticalPoint=1/(2*log(2)*alpha)-1./G;
    
    %compute powers
    for k=1:length(G)
        p(k)=min(lb(k),max(ub(k),criticalPoint(k)));
    end
    
    %compute F(alpha)=sumRate-alpha*powerConsumption
    R_P=(1/(2*length(G)))*sum(log2(1+G.*p'))-alpha*(sum(p)/length(G)+Pc);
    alpha=(1/(2*length(G)))*sum(log2(1+G.*p'))/(sum(p)/length(G)+Pc);
    

end

end