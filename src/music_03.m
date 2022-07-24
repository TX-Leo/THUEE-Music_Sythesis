clear;
close all;
clc;
for i = 1 : 7
    if (i ==1)
        tunes = zeros([7 4]);
        tunes(3, 1) = 220;                                          %f major-> 1 corresponding to f
        tunes(3, 2) = 440;
        tunes(3, 3) = 880;
        frequency_diff = 2^(1/12).^[-4, -2, 0, 2, 3, 5, 7]';        %The semitone frequency multiplier is 2^(1/12)
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
length_of_beat = 0.5;                                               %The length of time in a beat

bass = @(x) x;                                                      %low pitch
alto = @(x) x + 7;                                                  %mid pitch
treble = @(x) x + 14;                                               %high pitch
pause = @(x) 22;                                                    %pause

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

tone_sampling = [];
padding = 1;
last_padding = 0;
for i = 1 : 1 : length(song_length)
    f = tunes(song_pitch(i));                                                               %The pitch of each tone
    length_of_each_tone = song_length(i) * length_of_beat;                                  %The length of each tone
    length_of_each_tone_padding = frequency_sampling * length_of_each_tone * padding;
    t = linspace(0, length_of_each_tone * padding - 1 / frequency_sampling, length_of_each_tone_padding)';
    tone_sampling_temp = my_envelope(t/length_of_each_tone) .* sin(2 * pi * f * t);         %These sounds are represented by a sinusoidal signal with an amplitude of 1 and sampling frequency of 8kHz(after padding)
    if (last_padding == 0)
        tone_sampling = [tone_sampling; tone_sampling_temp];
    else
        tone_sampling = [tone_sampling(1:end-last_padding); tone_sampling(end-last_padding)+tone_sampling_temp(1:last_padding); tone_sampling_temp(last_padding+1:end)];
    end
    last_padding = round(length_of_each_tone_padding - frequency_sampling * length_of_each_tone);
end

% Raise an octave
subplot(3,2,1);
plot([0 : 500 - 1] / frequency_sampling*2, tone_sampling(1 : 500));
title("treble sound waveform");
subplot(3,2,2);
plot([0 : length(tone_sampling) - 1] / frequency_sampling * 2, tone_sampling);              %draw sound waveform
title("treble sound waveform");
% sound(tone_sampling, frequency_sampling * 2);                                               %play sampled tones
audiowrite('music_03_treble.wav', tone_sampling, frequency_sampling * 2);                   %store .wav

% Lower an octave
subplot(3,2,3);
plot([0 : 500 - 1] / frequency_sampling/2, tone_sampling(1 : 500));
title("bass sound waveform");
subplot(3,2,4);
plot([0 : length(tone_sampling) - 1] / frequency_sampling / 2, tone_sampling);              %draw sound waveform
title("bass sound waveform");
% sound(tone_sampling, frequency_sampling / 2);                                               %play sampled tones
audiowrite('music_03_bass.wav', tone_sampling, frequency_sampling / 2);                     %store .wav

%Raise half an octave
subplot(3,2,5);
tone_sampling_resample = resample(tone_sampling, round(frequency_sampling ./ 2.^ (1/12)), frequency_sampling);
plot([0 : 500 - 1] / frequency_sampling, tone_sampling_resample(1 : 500));
title("treble half sound waveform");
subplot(3,2,6);
plot([0 : length(tone_sampling_resample) - 1] / frequency_sampling, tone_sampling_resample);%draw sound waveform
title("treble half sound waveform");
sound(tone_sampling_resample, frequency_sampling);                                          %play sampled tones
audiowrite('music_03_treble_half.wav', tone_sampling_resample, frequency_sampling);         %store .wav

%envelope fuction
function [y, z] = my_envelope(x)
    z = 1;
    y = (x >= 0 & x < 0.1) .* ( (-1./(0.1.^2)).* (x-0.1).^2+ 1) + (x >= 0.1 & x < 1.5) .* ((((1 + (1.5-0.1)./(0.1-1))-1)./(1.5-0.1).^2) .* (x -0.1).^2 + 1) + (x >= 1.5 & x < 1) .* ((1.5-0.1)./(1.5-1) .*(((1 + (1.5-0.1)./(0.1-1))-1)./(1.5-0.1).^2)) .* (x -1).^2;
end