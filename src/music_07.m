clear;
close all;
clc;
load('resource/Guitar.MAT');
% Name             Size            Bytes  Class     Attributes
% realwave       243x1              1944  double
% wave2proc      243x1              1944  double

frequency_sampling = 8000;

subplot(3, 1, 1);
plot([0 : length(realwave)-1] / frequency_sampling, realwave);
title('realwave-Real guitar sound');
xlabel('t(8khz sampled)');
ylabel('A');

subplot(3, 1, 2);
plot([0 : length(wave2proc)-1] / frequency_sampling, wave2proc);
title('wave2proc-Guitar sounds to be dealt with');
xlabel('t(8khz sampled)');
ylabel('A');

subplot(3, 1, 3);
realwave_sampling_resample = my_resample(realwave);
    % realwave_resample = resample(realwave,100,1);
    % len = round(length(realwave_resample)/10);
    % realwave_sampling_resample = zeros([len,1]);
    % for i = 1:1:10
    %     realwave_sampling_resample = realwave_sampling_resample+realwave_resample((i-1)*len+1:i*len);
    % end

    % realwave_sampling_resample = [realwave_sampling_resample/10; realwave_sampling_resample/10];
    % realwave_sampling_resample = [realwave_sampling_resample; realwave_sampling_resample; realwave_sampling_resample; realwave_sampling_resample; realwave_sampling_resample];
    % realwave_sampling_resample = resample(realwave_sampling_resample, 1, 100);
plot([0 : length(realwave_sampling_resample)-1] / frequency_sampling, realwave_sampling_resample );
title('modifiedwave-modified guitar sound by realwave');
xlabel('t(8khz sampled)');
ylabel('A');
sound(realwave_sampling_resample, frequency_sampling);
audiowrite('music_07_modified_realwave.wav', realwave_sampling_resample, frequency_sampling);                           %store .wav

function sampling_resample = my_resample(x)
    resample_temp = resample(x,100,1);
    len = round(length(resample_temp)/10);
    sampling_resample = zeros([len,1]);
    for i = 1:1:10
        sampling_resample = sampling_resample+resample_temp((i-1)*len+1:i*len);
    end
    sampling_resample = [sampling_resample/10; sampling_resample/10];
    sampling_resample = [sampling_resample; sampling_resample; sampling_resample; sampling_resample; sampling_resample];
    sampling_resample = resample(sampling_resample, 1, 100);
end
