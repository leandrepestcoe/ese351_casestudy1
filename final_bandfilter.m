function [filter] = final_bandfilter(x,t_new)

% initial bandpass filter 
R_high0 = 1000;
C_high0 = 1.59e-4; % 1 Hz
R_low0 = 1000;
C_low0 = 7.96e-7; %200 Hz
y_high0 = lsim([1 0],[1 1/(R_high0*C_high0)], x, t_new);
y_low0 = lsim(1/(R_low0*C_low0),[1 1/(R_low0*C_low0)],x,t_new);
y_band0 = y_low0+y_high0;

% bandpass filter 1
R_high1 = 1000;
C_high1 = 7.96e-7; %200 Hz
R_low1 = 1000;
C_low1 = 2.89e-7; %550 Hz
y_low1 = lsim(1/(R_low1*C_low1),[1 1/(R_low1*C_low1)],x,t_new);
y_high1 = lsim([1 0],[1 1/(R_high1*C_high1)], x, t_new);
y_band1 = y_low1+y_high1;

% bandpass filter 2
R_high2 = 1000;
C_high2 = 2.89e-7; %550 Hz
R_low2 = 1000;
C_low2 = 1.77e-7; %900 Hz
y_low2 = lsim(1/(R_low2*C_low2),[1 1/(R_low2*C_low2)],x,t_new);
y_high2 = lsim([1 0],[1 1/(R_high2*C_high2)], x, t_new);
y_band2 = y_low2+y_high2;

% bandpass filter 3
R_high3 = 1000;
C_high3 = 1.77e-7; %900 Hz
R_low3 = 1000;
C_low3 = 1.27e-7; %1250 Hz
y_low3 = lsim(1/(R_low3*C_low3),[1 1/(R_low3*C_low3)],x,t_new);
y_high3 = lsim([1 0],[1 1/(R_high3*C_high3)], x, t_new);
y_band3 = y_low3+y_high3;

% bandpass filter 4
R_high4 = 1000;
C_high4 = 1.27e-7; %1250 Hz
R_low4 = 1000;
C_low4 = 9.95e-8; %1600 Hz
y_low4 = lsim(1/(R_low4*C_low4),[1 1/(R_low4*C_low4)],x,t_new);
y_high4 = lsim([1 0],[1 1/(R_high4*C_high4)], x, t_new);
y_band4 = y_low4+y_high4;

% bandpass filter 5
R_high5 = 1000;
C_high5 = 9.95e-8; %1600 Hz
R_low5 = 1000;
C_low5 = 2.65e-8; %6000 Hz
y_low5 = lsim(1/(R_low5*C_low5),[1 1/(R_low5*C_low5)],x,t_new);
y_high5 = lsim([1 0],[1 1/(R_high5*C_high5)], x, t_new);
y_band5 = y_low5+y_high5;


% gains for freq bands
gain0 = 0.01; %1-200 Hz
gain1 = 0.01; %200-550 Hz
gain2 = 0.01; %550-900 Hz
gain3 = 0.5; %900-1250 Hz
gain4 = 1; %1250-1600 Hz
gain5 = .3; %1600-6000 Hz

filter = (gain0*y_band0)+(gain1*y_band1)+(gain2*y_band2)+(gain3*y_band3)+(gain4*y_band4)+(gain5*y_band5); %six bandpass filters in parallel
%filter = y_band1;
end

