function output = func_Q_c(clk_bias,ISB_E,ISB_B,clk_drift)
% INPUT: standard deviation of GPS clock bias, 
%        Inter-System clock bias of Galilio, BeiDou compared to GPS clock
%        and GPS clock drift
    output = diag([clk_bias^2; ISB_E^2; ISB_B^2; clk_drift^2]);
end