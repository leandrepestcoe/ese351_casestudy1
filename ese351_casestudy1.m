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

f = 100;
w = 2*pi*f;
delta_t = 1/44100;
t_new = (0:delta_t:0.2);
%x = exp(j*w*t_new);
x = xv; %input is audio data...
R_low = 1000;
C_low = 0.00002;
R_high = 1000;
C_high = 0.00002;

% calculate output for low pass filter
y_low = lsim(1/(R_low*C_low),[1 1/(R_low*C_low)],x,t_new);
figure
hold on;
plot(t_new,x);
plot(t_new,y_low);
legend('input','output');
hold off;

% calculate output for highpass filter 
y_high = lsim([1 0],[1 1/(R_high*C_high)], x, t_new);
figure
hold on;
plot(t_new,x);
plot(t_new,y_high);
legend('input','output');
hold off