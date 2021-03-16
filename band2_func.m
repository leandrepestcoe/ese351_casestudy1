function [y_band2] = band2_func(x,t_new)
% bandpass filter 2
R_high2 = 1000;
C_high2 = 2.89e-7; %550 Hz
R_low2 = 1000;
C_low2 = 1.77e-7; %900 Hz
y_low2 = lsim(1/(R_low2*C_low2),[1 1/(R_low2*C_low2)],x,t_new);
y_high2 = lsim([1 0],[1 1/(R_high2*C_high2)], x, t_new);
y_band2 = y_low2+y_high2;
end

