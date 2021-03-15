%% Case Study 1
% Leandre Pestcoe and Julianne Wegmann

%% Load 'noisyviolin.mat' and Play Recording

[xv,xvfs] = audioread('violindirty.wav');
fs = xvfs; 
%sound(xv,fs)

%% Compute and Plot fft

f = [0:length(xv)/2]*fs/length(xv);
XV = fft(x);
P2 = abs(XV/length(xv));
P1 = P2(1:length(xv)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure, plot(f,P1);
xlabel('f, Hz')
ylabel('|X(f)|')

%% Define Bandpass Filter with 6 Freq Bands

delta_t = 0.1;
t_new = (0:delta_t:(length(xv)-1)*delta_t)'; %define time vector
x = xv; %set input to audio signal

final_filter = final_bandfilter(x,t_new);

sound(final_filter,fs);

%freq_response = freqs(final_filter,x,round(2*pi*60));
%plot(freq_response);
%freqz(final_filter);

%% Compute/Plot Freq Response

T = 0.002;
fs = 44100;
delta_t = 1/fs;
f_range = logspace(1,log10(20000),600);
t = (0:delta_t:10*T);
H = zeros(size(f_range));
x = exp(j*2*pi*100*t);

for i = 1:length(f_range)
    f = f_range(i);
    x = exp(j*2*pi*f*t);
    y = final_bandfilter(x,t);
    H(i) = y(end)/x(end); %use end because we want to analyze the steady state part 
end

gain_mag = 20*log10(abs(H)); %in dB units
gain_phase = (angle(H)/pi);

figure();
sgtitle('Bode Plot for Frequency Response of All Bands');
subplot(2,1,1);
semilogx(f_range, gain_mag);
xlim([0 10000]);
title('Magnitude of Gain');
xlabel('Frequency'); ylabel('dB');
subplot(2,1,2);
semilogx(f_range, gain_phase);
title('Phase of Gain');
xlabel('Frequency'); ylabel('Radians');
xlim([0 10000]);

%% Multi-band equalizer using bandpass func

band0 = bandpass(x, [50 200], fs);
band1 = bandpass(x, [200 550], fs);
band2 = bandpass(x, [550 900], fs);
band3 = bandpass(x, [900 1250], fs);
band4 = bandpass(x, [1250 1600], fs);
band5 = bandpass(x, [1600 2000], fs);
final_band = 0.1*band0+0.1*band1+2*band2+2*band3+band4+0.1*band5;
sound(final_band,fs);

%figure();
%pspectrum(final_band,fs);

%% frequency spectrum

y_out = fft(final_band);
L = length(final_band);

P2_out = abs(y_out/L);
P1_out = P2_out(1:L/2+1);
P1_out(2:end-1) = 2*P1_out(2:end-1);

f = fs*(0:(L/2))/L;

figure();
plot(f,P1_out) 
title('Single-Sided Amplitude Spectrum of x_out(t)')
xlabel('f (Hz)')
ylabel('|y_out(f)|')

%% experimenting with filter designer

% load the filter from the .mat file into the workspace
%good_filter = load('60 Hz Bandstop Filter.mat');

% apply the filter to the audio signal
%x_out = filter(good_filter.mat,x); % SOUNDS GREAT
x_out1 = filter(bandstop60,x);
x_out2 = filter(bandstop1900,x_out1);
x_out3 = filter(bandstop_low,x); % cuts out the low humming 
x_out4= filter (bandstop200_300,x);
finalfilter = 0.5*x_out3;
%x_out=x_out1+x_out2;
sound(finalfilter,xvfs);

% plot filtered signal
% figure();
% plot(t,x_out);
% %xlim([0 0.1]);
% title('filtered audio signal');
% xlabel('t');
% ylabel('x(t)');

%% 
% bode plot for magnitude vs. frequency
freq_range = logspace(1,log10(20000),10);
gains = zeros(size(freq_range));
T = 0.002;
dt = 20*T;
fs = 44100;
t = 0:(1/fs):dt;

for k = 1:length(freq_range);
    freq = freq_range(k);
    comp_exp = exp(j*freq*2*pi*t);
    band0 = bandpass(comp_exp, [50 200], fs);
    band1 = bandpass(comp_exp, [200 550], fs);
    band2 = bandpass(comp_exp, [550 900], fs);
    band3 = bandpass(comp_exp, [900 1250], fs);
    band4 = bandpass(comp_exp, [1250 1600], fs);
    band5 = bandpass(comp_exp, [1600 2000], fs);
    final_band = 0.1*band0+0.1*band1+band2+band3+band4+0.01*band5;
    filter_gain = final_band(end)/comp_exp(end);
    gains(k) = filter_gain;
end
% for k = 1:length(freq_range);
%     R_low0 = 1000;
%     C_low0 = 1.59e-4; %1 Hz
%     R_high0 = 1000;
%     C_high0 = 2.65e-6; % 60 Hz
%     y_high0 = lsim([1 0],[1 1/(R_high0*C_high0)], x, t_new);
%     y_low0 = lsim(1/(R_low0*C_low0),[1 1/(R_low0*C_low0)],x,t_new);
%     y_band0 = y_low0+y_high0;
%      final_band = y_band0;
%      filter_gain = final_band(end)/comp_exp(end);
%      gains(k) = filter_gain;
% end 

dB = 20*log10(abs(gains));
figure()
subplot(2,1,1);
semilogx(freq_range,dB);
title('Magnitude of H(w) for final filter');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;