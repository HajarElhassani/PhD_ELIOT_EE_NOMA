%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouput: G is a Tx2 matrix, each row t contains the channel gains for
%        link BS-user1 and BS-user2 at time t
%
% Inputs: T - time horizon for one run
%         sigma1, sigma2, var_h1, var_h2 - noise and channel link variances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [G] = channel_gain(T,sigma1,sigma2,var_h1,var_h2)

    
    % Wireless multi-path channels: Rayleigh distributed r.v. 
    
    h1=abs(sqrt(var_h1)*randn(T,1) + j*sqrt(var_h1)*randn(T,1));
    h2=abs(sqrt(var_h2)*randn(T,1) + j*sqrt(var_h2)*randn(T,1));
    
    % Normalized channel gains to the channel noise variance at the
    % receivers
    
    G=[h1.^2/sigma1, h2.^2/sigma2];
    
end