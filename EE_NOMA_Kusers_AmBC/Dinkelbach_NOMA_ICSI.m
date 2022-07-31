%------ Dinkelbach's algorithm -------% 

% input:          H -> channel gain of the direct link
%                 G -> channel gain of the backscattered+direct link
%                 H_error, G_error -> the channel gain with error
%                 A -V A=2^(2*Rmin) the SNR threshold
%                 Pmax -> maximum power budget 
%                 Pc -> circuit power (Pc)
%                 q -> the probability of backscattering (B=1)
% ouput:          alpha -> global energy efficiency (alpha)
%                 outage -> corresponding outage

function [alpha,outage] = Dinkelbach_NOMA_ICSI(q,H,G,H_error,G_error,A,Pmax,Pmin,Pc)

%---- initialization ----%
epsilon=1e-7;
alpha=1e-5;
R_P=1;
r=zeros(1,length(H));
idx=[];

% box constrainst : lb -> lower bound   ;  ub -> upper bound
lb=(A(1)-1)/H_error(1);
ub=1/prod(A(2:length(H_error)))*(Pmax-Pmin+lb*prod(A(2:length(H_error))));

%----- Dinkelbach's algorithm ----%
while(R_P>epsilon)
    %compute the optimal power p1
    p1=compute_EE_NOMA(q,alpha,A,H_error,G_error,Pc,ub,lb);
    %compute the optimal power allocation (theta_k)
    theta(1)=p1;
    theta(2)=(A(2)-1)/H_error(2)+theta(1)*A(2);
    
    if length(H_error)>2
        for k=3:length(H_error)
            B=0;
            for i=2:k-1
                B=B+(A(i)-1)/H_error(i)*prod(A(i+1:k));
            end
            theta(k)=(A(k)-1)/H_error(k)+B+theta(1)*prod(A(2:k)); %% \theta(k)=\sum_{i=1}^k p_i
        end
    end
    
    %------compute (gamma_{k \rightarrow i|0},gamma_{k \rightarrow i|1}) => direct(|0) and reflected (|1) links  with real channels-----%
    
    %compute(gamma_{k \rightarrow k|0},gamma_{k \rightarrow k|1})
    gamma0=zeros(length(H)); %to store gamma_{k \rightarrow k|0} in a matrix
    gamma1=zeros(length(H));  %to store gamma_{k \rightarrow k|1} in a matrix
    pass_afr=ones(length(H)-1,length(H)); %to store whether the conditions are met
    pass_bfr=ones(length(H)-1,length(H)); %to compare with pass_aft (if equal => constraints are met)
    
    gamma0(1,1)=(1+H(1)*theta(1)); % for user 1 (NOT decoded by any user)
    gamma1(1,1)=(1+G(1)*theta(1));
    
    for k=2:length(H)
        for i=1:k
            gamma0(i,k)=(1+H(i)*theta(k))/(1+H(i)*theta(k-1));
            gamma1(i,k)=(1+G(i)*theta(k))/(1+G(i)*theta(k-1));
        end 
    end
    
   
    % checking the constraints for SIC and QoS
    if (gamma0(1,1)>=A(1)-epsilon  && gamma1(1,1)>=A(1)-epsilon )
        pass_afr(1,1)=1;
    else
        pass_afr(1,1)=0;
    end
    for k=2:length(H)
        for i=1:k-1
            if (gamma0(i,k)>=gamma0(k,k)-epsilon && gamma1(i,k)>=gamma1(k,k)-epsilon && gamma0(k,k)>=A(k)-epsilon && gamma1(k,k)>=A(k)-epsilon)
                pass_afr(i,k)=1;
            else
                pass_afr(i,k)=0;
            end
        end
    end
    if (pass_afr==pass_bfr)
        r(1)=q/2*log2(gamma1(1,1))+(1-q)/2*log2(gamma0(1,1));
        for k=2:length(H)
            r(k)=q/2*log2(gamma1(k,k))+(1-q)/2*log2(gamma0(k,k));
        end
    else
        r=zeros(1,length(H));
    end
    
    %%%%%%%%%%%% check QoS constraints, we take A(end) since A is the same for all users %%%%%%%%%%%%%%%%
    Rmin=0.5*log2(A(end));
    idx=find(r<(Rmin-1e-5));
    
    %if at least one of the users doesn't satisfy his QoS constraints => outage
    if (isempty(idx))
        outage=0;
    else
        outage=1;
    end
    
    % compute the energy efficiency as the sum rate vs power consumption tradeoff
    R_P=sum(r)-alpha*(theta(end)+Pc); 
    % update alpha (energy efficiency as a ratio)
    alpha=sum(r)/(theta(end)+Pc);

end

end
