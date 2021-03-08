%% Case Study 1
% Leandre Pestcoe and Julianne Wegmann

%% load noisy violin data, then play recording
[xv,xvfs] = audioread('violindirty.wav');
fs = xvfs; 
sound(xv,fs)


%% compute and plot fft

f = [0:length(xv)-1]*fs/length(xv);
XV = fft(xv);
figure, plot(f,abs(XV));
xlabel('f, Hz')
ylabel('|X(f)|')

%%

f = 10;
w = 2*pi*f;
t_new = (0:delta_t:0.2);
x = exp(j*w*t_new);

%calculate output for low pass filter
y_low = lsim(1/(R*C),[1 1/(R*C)],x,t_new);
