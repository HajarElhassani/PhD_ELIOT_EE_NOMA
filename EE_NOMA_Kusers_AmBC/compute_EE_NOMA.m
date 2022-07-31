%% This function compute the optimal p1
% input      q -> the probability of backscattering (B=1)
%            alpha -> the tradeoff parameter
%            A -> the SNR threshold
%            H -> the channel gain of direct link
%            G -> the channel gain of the direct+backscattered link
%            Pc -> circuit power
%            ub -> the upper bound on p1
%            lb -> the lower bound on p1
% output     p1 -> the optimal power p1

function p1 = compute_EE_NOMA(q,alpha,A,H,G,Pc,ub,lb)

%-----------------------------------------------------------------%
% \theta_k=p_1*a_k+b_k                                                                %
%-----------------------------------------------------------------%

%number of users
K=length(H);
%initialize the vector theta
bk=zeros(K,1);
ak=zeros(K,1);

bk(1)=(A(1)-1)/H(1);
bk(2)=(A(2)-1)/H(2);
ak(2)=A(2);

%compute theta_k for k in {3,...,K}
for k=3:length(H)
    B=0;
    for i=2:(k-1)
        B=B+(A(i)-1)/H(i)*prod(A(i+1:k));   
    end
    bk(k)=B+(A(k)-1)/H(k); 
    ak(k)=prod(A(2:k));
end

%find the optimal power p1 using fmincon 
fun=@(p1) -(q/2*(sum(log2(1+G.*(p1*ak+bk)))-sum(log2(1+G(2:K).*(p1*ak(1:K-1)+bk(1:K-1)))))+(1-q)/2*log2(1+H(1)*p1)+sum((1-q)/2*log2(A(2:K)))-alpha*(p1*ak(K)+bk(K)+Pc));
Afun = [];
b = [];
Aeq = [];
beq = [];
x0=lb;
[p1,EE] = fmincon(fun,x0,Afun,b,Aeq,beq,lb,ub);

end