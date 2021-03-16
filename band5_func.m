function [y_band5] = band5_func(x,t_new)
% bandpass filter 5
R_high5 = 1000;
C_high5 = 9.95e-8; %1600 Hz
R_low5 = 1000;
C_low5 = 2.65e-8; %6000 Hz
y_low5 = lsim(1/(R_low5*C_low5),[1 1/(R_low5*C_low5)],x,t_new);
y_high5 = lsim([1 0],[1 1/(R_high5*C_high5)], x, t_new);
y_band5 = y_low5+y_high5;
end

