clear;
close all;
clc;
load('resource/Guitar.MAT');
% Name             Size            Bytes  Class     Attributes
% realwave       243x1              1944  double
% wave2proc      243x1              1944  double

frequency_sampling = 8000;

% wave2proc - one period of FFT
subplot(3, 1, 1);
fft1 = fft(wave2proc(1:round(end/10)));
plot([0 : length(fft1) - 1] / (round(length(wave2proc)/10)/frequency_sampling), abs(fft1));
title('wave2proc - one period of FFT');
xlabel('frequency');
ylabel('|w|');

% wave2proc - whole signal of FFT(10 periods)
subplot(3, 1, 2);
fft2 = fft(wave2proc);
plot([0 : length(fft2) - 1] / (length(wave2proc)/frequency_sampling), abs(fft2));
title('wave2proc - whole signal of FFT(10 periods)');
xlabel('frequency');
ylabel('|w|');

%wave2proc - 10*whole signal of FFT(100 periods
subplot(3, 1, 3);
for i = 1 : 1 : 10
    wave2proc = [wave2proc; wave2proc];
end
fft3 = fft(wave2proc);
plot([0 : length(fft3) - 1] / (length(wave2proc)/frequency_sampling), abs(fft3));
title('wave2proc - 10*whole signal of FFT(100 periods)');
xlabel('frequency');
ylabel('|w|');


%The intensity of each harmonic component
harmonic= my_get_harmonic(fft3);
disp([[0 : length(harmonic) - 1]', harmonic]);

function harmonic = my_get_harmonic(fft)
    delta = 10 * 2^10;
    n = [0 : length(fft) - 1];
    harmonic = abs(fft(mod(n, delta) == 0));
    harmonic = harmonic / harmonic(2);
end