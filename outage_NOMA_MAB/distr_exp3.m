
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: PhD student Hajar El Hassani, Anne Savard, E. Veronica Belmega
% ETIS Lab, UMR 8051, CY Cergy Paris Université, ENSEA, CNRS, F-95000,
% France
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the EXP3 probability distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouput: pr -> ax1 vector of the discrete probability distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%         weights -> ax1 vector of mutliplicative weights of EXP3
%         gamma -> coefficient between pure exploration with uniform
%         probability distribution and multiplicative weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pr = distr_exp3(weights,gamma)

    % normalization coefficient 
    weights_sum=sum(weights);
    
    % convex comination between uniform distribution and exponential
    % weights
    pr=(1-gamma)*(weights/weights_sum) + gamma/length(weights);

end

