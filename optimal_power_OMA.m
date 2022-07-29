%% This function computes the optimal power allocation for OMA and the corresponding sum rate and total power consumption
% input: A: 2^Rmin, G: channel gain, alpha: tradeoff parameter,
                    % Pmax: power budget, Pc: circuit power
%output: sum_Rate=sum rate, sum_Power=total power consumption

function [sum_Rate,sum_Power] = optimal_power_OMA(A,G,alpha,Pmax,Pc)
%lower bound of the feasible powers
p1_lb=(A(1)^3-1)/G(1);
p2_lb=(A(2)^3-1)/G(2);
p3_lb=(A(3)^3-1)/G(3);

%critical point
p1_c=1/(log(2)*alpha)-1/G(1);
p2_c=1/(log(2)*alpha)-1/G(2);
p3_c=1/(log(2)*alpha)-1/G(3);
% upper bound of the feasible powers
p1_ub=Pmax;
p2_ub=Pmax;
p3_ub=Pmax;

%the power allocation vector p=(p1,p2,p3)
p1=min(p1_ub,max(p1_lb,p1_c));
p2=min(p2_ub,max(p2_lb,p2_c));
p3=min(p3_ub,max(p3_lb,p3_c));

%compute the sum rate and power consumption
sum_Rate=(log2(1+G(1)*p1)+log2(1+G(2)*p2)+log2(1+G(3)*p3))/3;
sum_Power=(p1+p2+p3)/3+Pc;

end

