%% Case Study 1
% Leandre Pestcoe and Julianne Wegmann

%% load noisy violin data, then play recording
[xv,xvfs] = audioread('violindirty.wav');
fs = xvfs; 
sound(xv,fs)

%% compute and plot fft

f = [0:length(xv)/2]*fs/length(xv);
XV = fft(final_filter);
P2 = abs(XV/length(xv));
P1 = P2(1:length(xv)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure, plot(f,P1);
xlabel('f, Hz')
ylabel('|X(f)|')

%%
delta_t = 0.1;
t_new = (0:delta_t:(length(xv)-1)*delta_t)';
x = xv; %input is audio data...

% bandpass filter 1
%R_low1 = 1000;
%C_low1 = 0.00002;
R_high1 = 1000;
C_high1 = 7.96e-7;
%y_low1 = lsim(1/(R_low1*C_low1),[1 1/(R_low1*C_low1)],x,t_new);
y_high1 = lsim([1 0],[1 1/(R_high1*C_high1)], x, t_new);
%y_band1 = y_low1+y_high1;

% bandpass filter 2
R_low2 = 1000;
C_low2 = 7.96e-7; % 200 Hz
R_high2 = 1000;
C_high2 = 9.95e-8; % 1600 Hz
y_low2 = lsim(1/(R_low2*C_low2),[1 1/(R_low2*C_low2)],x,t_new);
y_high2 = lsim([1 0],[1 1/(R_high2*C_high2)], x, t_new);
y_band2 = y_low2+y_high2;

% bandpass filter 3
R_low3 = 1000;
C_low3 = 9.95e-8;
%R_high3 = 1000;
%C_high3 = 0.00002;
y_low3 = lsim(1/(R_low3*C_low3),[1 1/(R_low3*C_low3)],x,t_new);
%y_high3 = lsim([1 0],[1 1/(R_high3*C_high3)], x, t_new);
%y_band3 = y_low3+y_high3;

% final filter
final_filter = (0.01)*y_high1+y_band2+(0.01)*y_low3;
sound(final_filter,fs);

% % bandpass filter 4
% 
% R_low4 = 1000;
% C_low4 = 0.00002;
% R_high4 = 1000;
% C_high4 = 0.00002;
% y_low4 = lsim(1/(R_low4*C_low4),[1 1/(R_low4*C_low4)],x,t_new);
% y_high4 = lsim([1 0],[1 1/(R_high4*C_high4)], x, t_new);
% y_band4 = y_low4+y_high4;
% 
% % bandpass filter 5
% 
% R_low5 = 1000;
% C_low5 = 0.00002;
% R_high5 = 1000;
% C_high5 = 0.00002;
% y_low5 = lsim(1/(R_low5*C_low5),[1 1/(R_low5*C_low5)],x,t_new);
% y_high5 = lsim([1 0],[1 1/(R_high5*C_high5)], x, t_new);
% y_band5 = y_low5+y_high5;

