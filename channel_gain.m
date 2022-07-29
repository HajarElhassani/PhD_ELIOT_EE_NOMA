%% generate channel gains following the complex gaussian distribution with zero mean and variance 'ch_var'

function [G] = channel_gain(K,ch_var)
%generate channel
h=sqrt(ch_var)*randn(K,1) + j*sqrt(ch_var)*randn(K,1); 
%compute the channel gain
h2=abs(h).^2;
% sort the channel gains in descending (from the strong user to the weak user)
G=sort(h2,'descend');
end