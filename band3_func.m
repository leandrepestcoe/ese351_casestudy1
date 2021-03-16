function [y_band3] = band3_func(x,t_new)
% bandpass filter 3
R_high3 = 1000;
C_high3 = 1.77e-7; %900 Hz
R_low3 = 1000;
C_low3 = 1.27e-7; %1250 Hz
y_low3 = lsim(1/(R_low3*C_low3),[1 1/(R_low3*C_low3)],x,t_new);
y_high3 = lsim([1 0],[1 1/(R_high3*C_high3)], x, t_new);
y_band3 = y_low3+y_high3;
end

