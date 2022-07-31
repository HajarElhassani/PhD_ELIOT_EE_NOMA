%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ouput: threshold_OMA, the SNR threshold not the rate threshold
% Input: threshold_NOMA, 
%
% Since the rate thresholds must be equal between OMA and NOMA,
% then the SNR threshold_OMA is a function of the SNR threshold_NOMA and
% they are not identical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [threshold_OMA] = threshold_conversion_OMA(threshold_NOMA)

threshold_OMA = (1+threshold_NOMA)^2-1;

end

