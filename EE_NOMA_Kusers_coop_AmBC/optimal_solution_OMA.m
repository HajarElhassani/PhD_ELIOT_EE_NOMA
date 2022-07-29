%% This function compute the GEE according to Dinkelbach's algorithm

%% input: channel gain (G), A=2^(2*Rmin), maximum power budget (Pmax), circuit power (Pc)
%% ouput: global energy efficiency (alpha)
function alpha = optimal_solution_OMA(G,A,Pmax,Pc)

%% initialization parameters
epsilon=1e-7;
alpha=1e-5;
R_P=1;

%% lower and upper bound of feasible p1
lb=(A.^length(G)-1)./G;
ub=Pmax*ones(length(G),1);

%% Dinkelbach's algorithm
while(R_P>epsilon)   
    c=1/(2*log(2)*alpha)-1./G; % critical point
    
    %% compute the powers
    for k=1:length(G)
        p(k)=min(ub(k),max(lb(k),c(k)));
    end
    %% compute the function sum rate vs total power tradeoff
    R_P=1/(2*length(G))*sum(log2(1+G.*p'))-alpha*(sum(p)/length(G)+Pc);
    
    %% update alpha
    alpha=1/(2*length(G))*sum(log2(1+G.*p'))/(sum(p)/length(G)+Pc);
end
end