%% This function converts values from dbm to Watt

function [P_watt] = dbm_to_Watt(P_dbm)
P_watt=10.^((P_dbm-30)/10);
end

