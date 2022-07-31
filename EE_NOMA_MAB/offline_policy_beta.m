function mu_best = offline_policy_beta(users,beta,Pmax,Pc,threshold1,threshold2,sigma1,sigma2,var_h1,var_h2)

%calcul des rewards (gains) de chaque bras
for i=1:length(users)
mu(i)=expectedValue_reward_beta(users(i),beta,threshold1,threshold2,Pmax,sigma1,sigma2,var_h1,var_h2);
end

mu_a=(log2(1+threshold1)+log2(1+threshold2))*mu./(Pmax*beta+Pc);
%choisir le bras qui a la moyenne statique max du reward
[val1,idx1]=find(mu_a==max(mu_a));
m1=ceil(rand*length(idx1)); %Choose a max at random in case we have more than 2 arms
ind1=idx1(m1);

mu_best=mu_a(ind1);

end

