%% This function generate the channels between the BS and (users/backscatter)
function [H,I]=channelGain_BS(cor_user,alpha,sigma)
Dis_user_BS=abs(cor_user);
h=sqrt(Dis_user_BS.^(-alpha)./sigma);
[H,I]=sort(h,'descend'); % the SIC order
end

