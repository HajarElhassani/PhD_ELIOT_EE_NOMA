%% This function computes the optimal power allocation for NOMA and the corresponding sum rate and total power consumption
% input: A: 2^Rmin, G: channel gain, alpha: tradeoff parameter,
                    % Pmax: power budget, Pc: circuit power
%output: sum_Rate=sum rate, sum_Power=total power consumption

function [sum_Rate,sum_Power] = optimal_power_NOMA(A,G,alpha,Pmax,Pc)

p1_lb=(A(1)-1)/G(1); %lower bound of the feasible p1
p1_ub=(Pmax-A(3)*(A(2)-1)/G(2)-(A(3)-1)/G(3))/(A(2)*A(3)); %upper bound of the feasible p1
p1_c=1/(log(2)*alpha*A(2)*A(3))-1/G(1); % the critical point

% the power allocation vector p=(p1,p2,p3)
p1=max(p1_lb,min(p1_ub,p1_c)); % 
p2=(A(2)-1)*(p1+1/G(2));
p3=(A(3)-1)*(p1+p2+1/G(3));

%conmpute the rates
r1=log2(1+G(1)*p1);
r2=log2(1+G(2)*p2/(1+G(2)*p1));
r3=log2(1+G(3)*p3/(1+G(3)*(p1+p2)));

%compute the sum rate and total power consumption
sum_Rate=r1+r2+r3;
sum_Power=p1+p2+p3+Pc;

end

