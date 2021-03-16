%% Case Study 1
% Leandre Pestcoe and Julianne Wegmann

%% Load 'noisyviolin.mat' and Play Recording

[xv,xvfs] = audioread('MATLAB .wav');
fs = xvfs; 
sound(xv,fs)

%% Define Bandpass Filter with 6 Freq Bands

delta_t = 0.1;
t_new = (0:delta_t:(length(xv)-1)*delta_t)'; %define time vector
x = xv; %set input to audio signal

final_filter = final_bandfilter(x,t_new);

sound(final_filter,fs);

%% Define Filter using filterDesigner
x=xv;
%load the filters from the .mat files into the workspace
filter_obj0 = load('reject_60Hz.mat');
filter_60_hz = filter_obj0.reject_60Hz; %reject 60 Hz

filter_obj1 = load('bandpass1.mat');
bandpass1 = filter_obj1.bandpass1; %200-550Hz

filter_obj2 = load('bandpass2.mat');
bandpass2 = filter_obj2.bandpass2; %550-900 Hz

filter_obj3 = load('bandpass3.mat');
bandpass3 = filter_obj3.bandpass3; %900-1250 Hz

filter_obj4 = load('bandpass4.mat');
bandpass4 = filter_obj4.bandpass4; %1250-1600 Hz

%filter audio
x_out0 = filter(filter_60_hz,x); %doesn't need a gain b/c it just gets rid of 60Hz
x_out1 = filter(bandpass1,x);
x_out2 = filter(bandpass2,x);
x_out3 = filter(bandpass3,x);
x_out4 = filter(bandpass4,x);

gain1 = 0.01;
gain2 = 1;
gain3 = 1;
gain4 = 4;

x_out = x_out0+(x_out1*gain1)+(x_out2*gain2)+(x_out3*gain3)+(x_out4*gain4);
sound(x_out,fs);

%% Compute/Plot fft

f = [0:length(xv)/2]*fs/length(xv);
XV = fft(x_out);
P2 = abs(XV/length(xv));
P1 = P2(1:length(xv)/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure, plot(f,P1);
xlim([0 8000]);
title('FFT of Filtered ViolinDirty Audio (Filter Designer)');
xlabel('f, Hz')
ylabel('|X(f)|')


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
    y = filter(filter_60_hz,x);
    %y = filter(bandpass4,x);
    %y = final_bandfilter(x,t);
    H(i) = y(end)/x(end); %use end because we want to analyze the steady state part 
end

gain_mag0 = 20*log10(abs(H)); %in dB units
gain_phase0 = (angle(H)/pi);

for i = 1:length(f_range)
    f = f_range(i);
    x = exp(j*2*pi*f*t);
    y = filter(bandpass1,x);
    %y = filter(bandpass4,x);
    %y = final_bandfilter(x,t);
    H(i) = y(end)/x(end); %use end because we want to analyze the steady state part 
end

gain_mag1 = 20*log10(abs(H)); %in dB units
gain_phase1 = (angle(H)/pi);

for i = 1:length(f_range)
    f = f_range(i);
    x = exp(j*2*pi*f*t);
    y = filter(bandpass2,x);
    %y = filter(bandpass4,x);
    %y = final_bandfilter(x,t);
    H(i) = y(end)/x(end); %use end because we want to analyze the steady state part 
end

gain_mag2 = 20*log10(abs(H)); %in dB units
gain_phase2 = (angle(H)/pi);

for i = 1:length(f_range)
    f = f_range(i);
    x = exp(j*2*pi*f*t);
    y = filter(bandpass3,x);
    %y = filter(bandpass4,x);
    %y = final_bandfilter(x,t);
    H(i) = y(end)/x(end); %use end because we want to analyze the steady state part 
end

gain_mag3 = 20*log10(abs(H)); %in dB units
gain_phase3 = (angle(H)/pi);

for i = 1:length(f_range)
    f = f_range(i);
    x = exp(j*2*pi*f*t);
    %y = filter(filter_60_hz,x);
    y = filter(bandpass4,x);
    %y = final_bandfilter(x,t);
    H(i) = y(end)/x(end); %use end because we want to analyze the steady state part 
end

gain_mag4 = 20*log10(abs(H)); %in dB units
gain_phase4 = (angle(H)/pi);

for i = 1:length(f_range)
    f = f_range(i);
    x = exp(j*2*pi*f*t);
    y = filter(filter_60_hz,x)+filter(bandpass1,x)+filter(bandpass2,x)*1+filter(bandpass3,x)*1+filter(bandpass4,x);

    %y = filter(x_out,x);
    %y = filter(bandpass4,x);
    %y = final_bandfilter(x,t);
    H(i) = y(end)/x(end); %use end because we want to analyze the steady state part 
end

gain_mag = 20*log10(abs(H)); %in dB units
gain_phase = (angle(H)/pi);



figure();
sgtitle('Magnitude of Frequency Response for Individual Bands');
subplot(2,3,1);
semilogx(f_range, gain_mag0);
xlim([0 10000]);
title('Magnitude of Bandstop for 60Hz');
xlabel('Frequency (Hz)'); ylabel('dB');

subplot(2,3,2);
semilogx(f_range, gain_mag1);
xlim([0 10000]);
title('Magnitude of Bandpass 1');
xlabel('Frequency (Hz)'); ylabel('dB');

subplot(2,3,3);
semilogx(f_range, gain_mag2);
xlim([0 10000]);
title('Magnitude of Bandpass 2');
xlabel('Frequency (Hz)'); ylabel('dB');

subplot(2,3,4);
semilogx(f_range, gain_mag3);
xlim([0 10000]);
title('Magnitude of Bandpass 3');
xlabel('Frequency (Hz)'); ylabel('dB');

subplot(2,3,5);
semilogx(f_range, gain_mag4);
xlim([0 10000]);
title('Magnitude of Bandpass 4');
xlabel('Frequency (Hz)'); ylabel('dB');

figure 
sgtitle('Phase of Frequency Response for Individual Bands');
subplot(2,3,1);
semilogx(f_range, gain_phase0);
title('Phase of Bandstop for 60 Hz');
xlabel('Frequency (Hz)'); ylabel('Radians');
xlim([0 10000]);

subplot(2,3,2);
semilogx(f_range, gain_phase1);
title('Phase of Bandpass 1');
xlabel('Frequency (Hz)'); ylabel('Radians');
xlim([0 10000]);

subplot(2,3,3);
semilogx(f_range, gain_phase2);
title('Phase of Bandpass 2');
xlabel('Frequency (Hz)'); ylabel('Radians');
xlim([0 10000]);

subplot(2,3,4);
semilogx(f_range, gain_phase3);
title('Phase of Bandpass 3');
xlabel('Frequency (Hz)'); ylabel('Radians');
xlim([0 10000]);

subplot(2,3,5);
semilogx(f_range, gain_phase4);
title('Phase of Bandpass 4');
xlabel('Frequency (Hz)'); ylabel('Radians');
xlim([0 10000]);

figure 
sgtitle('Bode Plot of Frequency Response for Final Band (with no gain)');
subplot(2,1,1)
semilogx(f_range, gain_mag);
title('Magnitude');
xlabel('Frequency (Hz)'); ylabel('dB');
xlim([0 10000]);

subplot(2,1,2)
semilogx(f_range, gain_phase);
title('Phase');
xlabel('Frequency (Hz)'); ylabel('Radians');
xlim([0 10000]);
%% Compute/Plot Impulse Response

T = 0.002;
f = 44100;
delta_t = 1/f;
t = (0:delta_t:T);
impulse = zeros(1,length(t));
impulse(1)=1; %define impulse func

%compute impulse response
filter_im = final_bandfilter(impulse,t);  

figure();
plot(t, filter_im); %plot impulse response
title('Impulse Response of Final Filter');
xlabel('Time(s)');
