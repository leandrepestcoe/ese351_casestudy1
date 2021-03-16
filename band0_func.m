function [y_band0] = band0_func(x,t_new)
% initial bandpass filter 
R_high0 = 1000;
C_high0 = 1.59e-4; % 1 Hz
R_low0 = 1000;
C_low0 = 7.96e-7; %200 Hz
y_high0 = lsim([1 0],[1 1/(R_high0*C_high0)], x, t_new);
y_low0 = lsim(1/(R_low0*C_low0),[1 1/(R_low0*C_low0)],x,t_new);
y_band0 = y_low0+y_high0;
end

