function [arm,arm_ind] = argmax(arms,empirical_mean,delta)

% Upper confidence bound at time t for all arms
Q=empirical_mean+delta;
% find the best arm maximizing Q
[val,idx]=find(Q==max(Q));
m=ceil(rand*length(idx));
arm_ind=idx(m);
arm=arms(arm_ind,:);

end

