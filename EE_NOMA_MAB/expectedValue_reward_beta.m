function expectedValue = expectedValue_reward_beta(user,beta,threshold1,threshold2,Pmax,sigma1,sigma2,var_h1,var_h2)

if user==1 
    P=[0.25*Pmax*beta 0.75*Pmax*beta];
    value1=sigma1*threshold1/P(1);
    value2=sigma1*threshold2/(P(2)-threshold2*P(1));
    value3=sigma2*threshold2/(P(2)-threshold2*P(1));
    variance1=var_h1;
    variance2=var_h2;
else
    P=[0.75*Pmax*beta 0.25*Pmax*beta];
    value1=sigma2*threshold2/P(2);
    value2=sigma2*threshold1/(P(1)-threshold1*P(2));
    value3=sigma1*threshold1/(P(1)-threshold1*P(2));
    variance1=var_h2;
    variance2=var_h1;
end
if value2<0
    expectedValue=0;
else
    expectedValue=exp(-max(value1,value2)/(2*variance1))*exp(-value3/(2*variance2));
end

end 