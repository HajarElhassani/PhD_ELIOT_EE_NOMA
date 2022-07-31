%% This function compute the exploration term in UCB

function [delta] = upper_bound(alpha,t,nb_sel_arm)

delta=sqrt(alpha*log(t)./(2*nb_sel_arm));

end

