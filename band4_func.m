function [y_band4] = band4_func(x,t_new)
% bandpass filter 4
R_high4 = 1000;
C_high4 = 1.27e-7; %1250 Hz
R_low4 = 1000;
C_low4 = 9.95e-8; %1600 Hz
y_low4 = lsim(1/(R_low4*C_low4),[1 1/(R_low4*C_low4)],x,t_new);
y_high4 = lsim([1 0],[1 1/(R_high4*C_high4)], x, t_new);
y_band4 = y_low4+y_high4;
end

