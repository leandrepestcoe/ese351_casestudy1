%% Case Study 1
% Leandre Pestcoe and Julianne Wegmann

%% load noisy violin data, then play recording
[xv,xvfs] = audioread('violindirty.wav');
fs = xvfs; 
sound(xv,fs)

%% compute and plot fft

f = [0:length(xv)/2]*fs/length(xv);
XV = fft(final_band);
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

% initial bandpass filter 
R_low0 = 1000;
C_low0 = 1.59e-4; %1 Hz
R_high0 = 1000;
C_high0 = 7.96e-7; %200 Hz
y_high0 = lsim([1 0],[1 1/(R_high0*C_high0)], x, t_new);
y_low0 = lsim(1/(R_low0*C_low0),[1 1/(R_low0*C_low0)],x,t_new);
%y_band0 = lsim(1/(R_high0*C_high0),[1 1/(R_high0*C_high0)],x,t_new);
y_band0 = y_low0+y_high0;

% bandpass filter 1
R_low1 = 1000;
C_low1 = 7.96e-7; %200 Hz
R_high1 = 1000;
C_high1 = 2.89e-7; %550 Hz
y_low1 = lsim(1/(R_low1*C_low1),[1 1/(R_low1*C_low1)],x,t_new);
y_high1 = lsim([1 0],[1 1/(R_high1*C_high1)], x, t_new);
%y_band1 = lsim(1/(R_high1*C_high1),[1 1/(R_high1*C_high1)],x,t_new);
y_band1 = y_low1+y_high1;

% bandpass filter 2
R_low2 = 1000;
C_low2 = 2.89e-7; %550 Hz
R_high2 = 1000;
C_high2 = 1.77e-7; %900 Hz
y_low2 = lsim(1/(R_low2*C_low2),[1 1/(R_low2*C_low2)],x,t_new);
y_high2 = lsim([1 0],[1 1/(R_high2*C_high2)], x, t_new);
%y_band2 = lsim(1/(R_high2*C_high2),[1 1/(R_high2*C_high2)],x,t_new);
y_band2 = y_low2+y_high2;

% bandpass filter 3
R_low3 = 1000;
C_low3 = 1.77e-7; %900 Hz
R_high3 = 1000;
C_high3 = 1.27e-7; %1250 Hz
y_low3 = lsim(1/(R_low3*C_low3),[1 1/(R_low3*C_low3)],x,t_new);
y_high3 = lsim([1 0],[1 1/(R_high3*C_high3)], x, t_new);
%y_band3 = lsim(1/(R_high3*C_high3),[1 1/(R_high3*C_high3)],x,t_new);
y_band3 = y_low3+y_high3;

% bandpass filter 4

R_low4 = 1000;
C_low4 = 1.27e-7; %1250 Hz
R_high4 = 1000;
C_high4 = 9.95e-8; %1600 Hz
y_low4 = lsim(1/(R_low4*C_low4),[1 1/(R_low4*C_low4)],x,t_new);
y_high4 = lsim([1 0],[1 1/(R_high4*C_high4)], x, t_new);
%y_band4 = lsim(1/(R_high4*C_high4),[1 1/(R_high4*C_high4)],x,t_new);
y_band4 = y_low4+y_high4;

% bandpass filter 5

R_low5 = 1000;
C_low5 = 9.95e-8; %1600 Hz
R_high5 = 1000;
C_high5 = 7.96e-8; %2000 Hz
y_low5 = lsim(1/(R_low5*C_low5),[1 1/(R_low5*C_low5)],x,t_new);
y_high5 = lsim([1 0],[1 1/(R_high5*C_high5)], x, t_new);
%y_band5 = lsim(1/(R_high5*C_high5),[1 1/(R_high5*C_high5)],x,t_new);
y_band5 = y_low5+y_high5;


% final filter
gain0 = 0.1;
gain1 = 1;
gain2 = 1;
gain3 = 1;
gain4 = 1;
gain5 = 0.5;
final_filter = (gain0*y_band0)+(gain1*y_band1)+(gain2*y_band2)+(gain3*y_band3)+(gain4*y_band4)+(gain5*y_band5); %five bandpass filters in parallel
%sound(final_filter,fs);

%freq_response = freqs(final_filter,x,round(2*pi*60));
%plot(freq_response);
%freqz(final_filter);


%%
 band0 = bandstop(x,[1 200], fs);
band1 = bandpass(x,[200 550], fs);
band2 = bandpass(x, [550 900], fs);
band3 = bandpass(x, [900 1250], fs);
band4 = bandpass(x, [1250 1600], fs);
 band5 = bandstop(x, [1600 2000], fs);
final_band = .1*band0+0.9*band1+band2+band3+band4+.1*band5;
sound(final_band,fs);
figure
pspectrum(final_band,fs);


