%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the optimal policy by exhaustive search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouputs: mu_best -> the maximum energy efficiency (reward) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: %Pmax -> available transmit power at the base station
%         Pc -> circuit power
%         threshold 1 and 2 -> the SNR minimum thresholds
%         sigma1, sigma2, var_h1, var_h2 -> noise and channel links variances
%         of user 1 and user 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mu_best = optimal_NOMA(threshold1,threshold2,Pmax,Pc,sigma1,sigma2,var_h1,var_h2)

% the range of powers
P1=(1e-2:1e-2:Pmax);
P2=(1e-2:1e-2:Pmax);

for i=1:length(P1)
    for j=1:length(P2)
        if (P1(i)+P2(j)<=Pmax)
            if (P2(j)-threshold2*P1(i)>0)
                variance1=var_h1;
                variance2=var_h2;
                value1=sigma1*threshold1/P1(i);
                value2=sigma1*threshold2/(P2(j)-threshold2*P1(i));
                value3=sigma2*threshold2/(P2(j)-threshold2*P1(i));
                if value2<0
                    expectedValue(i,j)=0;
                else
                    expectedValue(i,j)=(log2(1+threshold1)+log2(1+threshold2))*exp(-max(value1,value2)/(2*variance1))*exp(-value3/(2*variance2))/(P1(i)+P2(j)+Pc);
                end
            end
            if (P1-threshold1*P2>0)
                variance1=var_h2;
                variance2=var_h1;
                value1=sigma2*threshold2/P2(j);
                value2=sigma2*threshold1/(P1(i)-threshold1*P2(j));
                value3=sigma1*threshold1/(P1(i)-threshold1*P2(j));
                if value2<0
                    expectedValue(i,j)=0;
                else
                    expectedValue(i,j)=(log2(1+threshold1)+log2(1+threshold2))*exp(-max(value1,value2)/(2*variance1))*exp(-value3/(2*variance2))/(P1(i)+P2(j)+Pc);
                end
            end
        else
            expectedValue(i,j)=0;
        end
    end   
end

mu_best=max(max(expectedValue));

end