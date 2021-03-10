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

%% define bandpass filter with 5 freq bands

delta_t = 0.1;
t_new = (0:delta_t:(length(xv)-1)*delta_t)'; %define time vector
x = xv; %set input to audio signal

% bandpass filter 1
R_low1 = 1000;
C_low1 = 7.96e-7; %200 Hz
R_high1 = 1000;
C_high1 = 3.32e-7; %480 Hz
y_low1 = lsim(1/(R_low1*C_low1),[1 1/(R_low1*C_low1)],x,t_new);
y_high1 = lsim([1 0],[1 1/(R_high1*C_high1)], x, t_new);
y_band1 = y_low1+y_high1;

% bandpass filter 2
R_low2 = 1000;
C_low2 = 3.32e-7; %480 Hz
R_high2 = 1000;
C_high2 = 2.09e-7; %760 Hz
y_low2 = lsim(1/(R_low2*C_low2),[1 1/(R_low2*C_low2)],x,t_new);
y_high2 = lsim([1 0],[1 1/(R_high2*C_high2)], x, t_new);
y_band2 = y_low2+y_high2;

% bandpass filter 3
R_low3 = 1000;
C_low3 = 2.09e-7; %760 Hz
R_high3 = 1000;
C_high3 = 1.53e-7; %1040 Hz
y_low3 = lsim(1/(R_low3*C_low3),[1 1/(R_low3*C_low3)],x,t_new);
y_high3 = lsim([1 0],[1 1/(R_high3*C_high3)], x, t_new);
y_band3 = y_low3+y_high3;

% bandpass filter 4

R_low4 = 1000;
C_low4 = 1.53e-7; %1040 Hz
R_high4 = 1000;
C_high4 = 1.21e-7; %1320 Hz
y_low4 = lsim(1/(R_low4*C_low4),[1 1/(R_low4*C_low4)],x,t_new);
y_high4 = lsim([1 0],[1 1/(R_high4*C_high4)], x, t_new);
y_band4 = y_low4+y_high4;

% bandpass filter 5

R_low5 = 1000;
C_low5 = 1.21e-7; %1320 Hz
R_high5 = 1000;
C_high5 = 9.95e-8; %1600 Hz
y_low5 = lsim(1/(R_low5*C_low5),[1 1/(R_low5*C_low5)],x,t_new);
y_high5 = lsim([1 0],[1 1/(R_high5*C_high5)], x, t_new);
y_band5 = y_low5+y_high5;


% final filter
gain1 = 1;
gain2 = 1;
gain3 = 1;
gain4 = 1;
gain5 = 1;
final_filter = (gain1*y_high1)+(gain2*y_band2)+(gain3*y_band3)+(gain4*y_band4)+(gain5*y_band5); %five bandpass filters in parallel
%sound(final_filter,fs);

freq_response = freqs(final_filter,x,round(2*pi*60));
plot(freq_response);

