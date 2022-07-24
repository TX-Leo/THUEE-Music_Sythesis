clear;
close all;
clc;
for i = 1 : 7
    if (i ==1)
        tunes = zeros([7 4]);
        tunes(3, 1) = 220;                                                  %f major-> 1 corresponding to f
        tunes(3, 2) = 440;
        tunes(3, 3) = 880;
        frequency_diff = 2^(1/12).^[-4, -2, 0, 2, 3, 5, 7]';           %The semitone frequency multiplier is 2^(1/12)
    end
    tunes(i, 1:end) = tunes(3, 1:end) .* frequency_diff(i);
end
% tunes=%[
% 175    349    698      0
% 196    392    784      0
% 220    440    880      0
% 247    494    988      0
% 262    523   1047      0
% 294    587   1175      0
% 330    659   1319      0
% ]

frequency_sampling = 8000;                                          %sampling frequency
length_of_beat = 0.5;

bass = @(x) x;                                                      %low pitch
alto = @(x) x + 7;                                                  %mid pitch
treble = @(x) x + 14;                                               %high pitch
pause = @(x) 22;

%The pitch of each tone
song_pitch =[...
treble(5); treble(5); treble(6);...                                 %5-  56
treble(2);...                                                       %2-  --
treble(1); treble(1); bass(6);...                                   %1-  16
treble(2)];
%The length of each tone
song_length =[...
1;0.5;0.5;...                                                       %5-  56
2;...                                                               %2-  --
1;0.5;0.5;...                                                       %1-  16
2];

%music_08;
%%%%%%%%%%%%%%%%%%%%music_08->to get variable harmonic%%%%%%%%%%%%%%%%%
load('resource/Guitar.MAT');
% Name             Size            Bytes  Class     Attributes
% realwave       243x1              1944  double
% wave2proc      243x1              1944  double

%wave2proc - 10*whole signal of FFT(100 periods
subplot(3, 1, 3);
for i = 1 : 1 : 10
    wave2proc = [wave2proc; wave2proc];
end
fft3 = fft(wave2proc);
% plot([0 : length(fft3) - 1] / (length(wave2proc)/frequency_sampling), abs(fft3));
% title('wave2proc - 10*whole signal of FFT(100 periods)');
% xlabel('frequency');
% ylabel('|w|');

%The intensity of each harmonic component
harmonic= my_get_harmonic(fft3);
disp([[0 : length(harmonic) - 1]', harmonic]);

tone_sampling = [];
padding = 1;
last_padding = 0;
for i = 1 : 1 : length(song_length)
    f = tunes(song_pitch(i));                                                               %The pitch of each tone
    length_of_each_tone = song_length(i) * length_of_beat;                                  %The length of each tone
    length_of_each_tone_padding = frequency_sampling * length_of_each_tone * padding;
    t = linspace(0, length_of_each_tone * padding - 1 / frequency_sampling, length_of_each_tone_padding)';
    tone_sampling_temp = zeros(size(t));
    for j = 1 : 1 : length(harmonic)
        tone_sampling_temp = tone_sampling_temp + harmonic(j) * sin(2 * pi * j * f * t);
    end
    tone_sampling_temp = my_envelope(t/length_of_each_tone) .* tone_sampling_temp;
    % tone_sampling_temp = my_envelope(t/length_of_each_tone) .* (sin(2 * pi * f * t) + 0.2 * sin(2 * pi * 2 * f * t) + 0.1 * sin(2 * pi * 3 * f * t)+0.4 * sin(2 * pi * 4 * f * t));         %These sounds are represented by a sinusoidal signal with an amplitude of 1 and sampling frequency of 8kHz(after padding)[The fundamental wave amplitude is 1, the second harmonic amplitude is 0:2, the third harmonic amplitude is 0:3]
    if (last_padding == 0)
        tone_sampling = [tone_sampling; tone_sampling_temp];
    else
        tone_sampling = [tone_sampling(1:end-last_padding); tone_sampling(end-last_padding)+tone_sampling_temp(1:last_padding); tone_sampling_temp(last_padding+1:end)];
    end
    last_padding = round(length_of_each_tone_padding - frequency_sampling * length_of_each_tone);
end

subplot(3,1,1);                                                     %draw envelope functiofn
fplot(@(t) my_envelope(t), [-padding+1, padding * 1.1]);
title("envelope function");
subplot(3,1,2);
plot([0 : 5000 - 1] / frequency_sampling, tone_sampling(1 : 5000));
title("sound waveform");
subplot(3,1,3);
plot([0 : length(tone_sampling) - 1] / frequency_sampling, tone_sampling);
title("sound waveform");

sound(tone_sampling, frequency_sampling);                       %play sampled tones
audiowrite('music_10.wav', tone_sampling, frequency_sampling);  %store .wav

function harmonic = my_get_harmonic(fft)
    delta = 10 * 2^10;
    n = [0 : length(fft) - 1];
    harmonic = abs(fft(mod(n, delta) == 0));
    harmonic = harmonic / harmonic(2);
end
function [y, z] = my_envelope(x)
    z = 1;
    y = (x >= 0 & x < 0.1) .* ( (-1./(0.1.^2)).* (x-0.1).^2+ 1) + (x >= 0.1 & x < 1.5) .* ((((1 + (1.5-0.1)./(0.1-1))-1)./(1.5-0.1).^2) .* (x -0.1).^2 + 1) + (x >= 1.5 & x < 1) .* ((1.5-0.1)./(1.5-1) .*(((1 + (1.5-0.1)./(0.1-1))-1)./(1.5-0.1).^2)) .* (x -1).^2;
end

