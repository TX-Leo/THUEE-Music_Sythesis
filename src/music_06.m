clear;
close all;
clc;
load('resource/Guitar.MAT');
% Name             Size            Bytes  Class     Attributes
% realwave       243x1              1944  double
% wave2proc      243x1              1944  double

frequency_sampling_guitar = 8000;

subplot(3, 1, 1);
t1 = linspace(0, length(realwave) / frequency_sampling_guitar - 1/frequency_sampling_guitar, length(realwave));
plot(t1, realwave);
title('realwave-Real guitar sound');
xlabel('t(8khz sampled)');
ylabel('A');
% sound(realwave, frequency_sampling_guitar);
audiowrite('music_06_realwave.wav', realwave, frequency_sampling_guitar);                           %store .wav

subplot(3, 1, 2);
t2 = linspace(0, length(wave2proc) / frequency_sampling_guitar - 1/frequency_sampling_guitar, length(wave2proc));
plot(t2, wave2proc);
title('wave2proc-Guitar sounds to be dealt with');
xlabel('t(8khz sampled)');
ylabel('A');
% sound(wave2proc, frequency_sampling_guitar);
audiowrite('music_06_wave2proce.wav', wave2proc, frequency_sampling_guitar);                        %store .wav

[tone_sampling, frequency_sampling] = audioread('resource/fmt.wav');
subplot(3, 1, 3);
t3 = linspace(0, length(tone_sampling) / frequency_sampling - 1/frequency_sampling, length(tone_sampling));
plot(t3, tone_sampling);
title('tone_sampling');
xlabel('t');
ylabel('A');
sound(tone_sampling, frequency_sampling);
audiowrite('music_06_fmt.wav', tone_sampling, frequency_sampling);                           %store .wav