%% This function compute the optimal power allocation vector using fmincon
% input      q -> the probability of backscattering (B=1)
%            alpha -> the tradeoff parameter
%            A -> the SNR threshold
%            H -> the channel gain of direct link
%            G -> the channel gain of the direct+backscattered link
%            Pc -> circuit power
% output     p -> the optimal power vector

function [p] = compute_EE_OMA(q,alpha,Pmax,A,H,G,Pc)

%upper bound and lower bound on the powers
lb=((A.^length(H)-1)./H);
ub=Pmax*ones(length(H),1);

fun=@(p) -(sum((q/(2*length(G)))*log2(1+G.*p')+((1-q)/(2*length(H)))*log2(1+H.*p'))-alpha*(sum(p)/length(H)+Pc));
A = [];
b = [];
Aeq = [];
beq = [];
x0=lb;
[p,EE] = fmincon(fun,lb',A,b,Aeq,beq,lb',ub');

end

