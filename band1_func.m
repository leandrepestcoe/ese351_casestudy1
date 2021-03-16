function [y_band1] = band1_func(x,t_new)
% bandpass filter 1
R_high1 = 1000;
C_high1 = 7.96e-7; %200 Hz
R_low1 = 1000;
C_low1 = 2.89e-7; %550 Hz
y_low1 = lsim(1/(R_low1*C_low1),[1 1/(R_low1*C_low1)],x,t_new);
y_high1 = lsim([1 0],[1 1/(R_high1*C_high1)], x, t_new);
y_band1 = y_low1+y_high1;
end

