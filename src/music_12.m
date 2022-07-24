function varargout = music_12(varargin)
% MUSIC_12 MATLAB code for music_12.fig
%      MUSIC_12, by itself, creates a new MUSIC_12 or raises the existing
%      singleton*.
%
%      H = MUSIC_12 returns the handle to a new MUSIC_12 or the handle to
%      the existing singleton*.
%
%      MUSIC_12('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MUSIC_12.M with the given input arguments.
%
%      MUSIC_12('Property','Value',...) creates a new MUSIC_12 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before music_12_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to music_12_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help music_12

% Last Modified by GUIDE v2.5 20-Jul-2022 22:44:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @music_12_OpeningFcn, ...
                   'gui_OutputFcn',  @music_12_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before music_12 is made visible.
function music_12_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to music_12 (see VARARGIN)

% Choose default command line output for music_12
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes music_12 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = music_12_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in add_audio.
function add_audio_Callback(hObject, eventdata, handles)
% hObject    handle to add_audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [file,path,~] = uigetfile('*.wav','请选择要识别的音频'); % Open the file selection dialog box
    if isequal(file,0)
       disp('User selected Cancel');
    else
       disp(['User selected ', fullfile(path,file)]);
    end
    file_path = fullfile(path,file);
    handles.audio_file_path = file_path;
    guidata(hObject, handles);

    output_string = file_path;
    set(handles.add_audio_file_path, 'String', output_string);


% --- Executes on button press in start_sampling.
function start_sampling_Callback(hObject, eventdata, handles)
% hObject    handle to start_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    file_path = handles.audio_file_path;
    % cmd = char("get_component(""" + file_path + """);");
    % disp(cmd);
    % [freq, comp] = evalin('base', cmd); % Computes MATLAB expressions in the specified workspace.
    [freq, comp] = get_component(file_path);
    handles.freq = freq;
    handles.comp = comp;
    guidata(hObject, handles);
    output_string = cell(length(freq), 1);
    for i = 1 : 1 : length(freq)
        output_string{i} = char(string(num2str(freq(i))) + ": [" + string(num2str(comp{i}')) + "]");
    end
    set(handles.harmonic_component, 'String', output_string);

% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in add_audio_file_path.
function add_audio_file_path_Callback(hObject, eventdata, handles)
% hObject    handle to add_audio_file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns add_audio_file_path contents as cell array
%        contents{get(hObject,'Value')} returns selected item from add_audio_file_path


% --- Executes during object creation, after setting all properties.
function add_audio_file_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to add_audio_file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in play_audio.
function play_audio_Callback(hObject, eventdata, handles)
% hObject    handle to play_audio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    file_path = handles.audio_file_path;
    [tone_sampling, frequency_sampling] = audioread(file_path);
    sound(tone_sampling, frequency_sampling);                           %play sampled tones
    axes(handles.audio_wave);
    plot([0 : length(tone_sampling) - 1] / frequency_sampling, tone_sampling);
    title("sound waveform");

% --- Executes on button press in compose_music.
function compose_music_Callback(hObject, eventdata, handles)
% hObject    handle to compose_music (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in import_music_table.
function import_music_table_Callback(hObject, eventdata, handles)
    % hObject    handle to import_music_table (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
% --- Executes on selection change in harmonic_component.
    [file,path,~] = uigetfile('*.xlsx','请选择要识别的表格'); % Open the file selection dialog box
    if isequal(file,0)
    disp('User selected Cancel');
    else
    disp(['User selected ', fullfile(path,file)]);
    end
    table_file_path = fullfile(path,file);
    output_string = table_file_path;
    set(handles.table_file_path, 'String', output_string);

    [song_pitch,song_length] = my_process_music_data_table(table_file_path);
    song = [song_pitch,song_length];
    set(handles.music_table,'Data',song);%将data中的文件以Data的形式设置在句柄为uitable1的表格中。

    handles.song_pitch = song_pitch;
    handles.song_length = song_length;
    guidata(hObject, handles);

% --- Executes on button press in input_parameter_and_generate_music.
function input_parameter_and_generate_music_Callback(hObject, eventdata, handles)
    % hObject    handle to input_parameter_and_generate_music (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    song_pitch = handles.song_pitch;
    song_length = handles.song_length;
    title_of_music = get(handles.title_of_music, 'String');
    main_tune = get(handles.tune, 'String');
    length_of_beat = str2num(get(handles.length_of_beat, 'String'));
    frequency_sampling = str2num(get(handles.frequency_sampling, 'String'));% frequency_sampling = 8000;
    harmonic = [1,0.2,0.1,0.4];
    tone_sampling = my_play_music(song_pitch,song_length, main_tune,length_of_beat,frequency_sampling,harmonic)

    handles.tone_sampling = tone_sampling;
    guidata(hObject, handles);

% --- Executes on button press in generate_wave_and_play_music.
function generate_wave_and_play_music_Callback(hObject, eventdata, handles)
    % hObject    handle to generate_wave_and_play_music (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    tone_sampling = handles.tone_sampling ;
    frequency_sampling = str2num(get(handles.frequency_sampling, 'String'));

    offset = str2num(get(handles.offset, 'String'));
    axes(handles.sound_wave);
    plot([0 : length(tone_sampling) - 1] / frequency_sampling * offset, tone_sampling);              %draw sound waveform
    title("sound waveform");
    sound(tone_sampling, frequency_sampling * offset);                                               %play sampled tones
    audiowrite('music_12_the_greatest_work.wav', tone_sampling, frequency_sampling * offset);                     %store .wav

function harmonic_component_Callback(hObject, eventdata, handles)
% hObject    handle to harmonic_component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns harmonic_component contents as cell array
%        contents{get(hObject,'Value')} returns selected item from harmonic_component


% --- Executes during object creation, after setting all properties.
function harmonic_component_CreateFcn(hObject, eventdata, handles)
% hObject    handle to harmonic_component (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function audio_wave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to audio_wave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate audio_wave

% --- Executes during object creation, after setting all properties.
function sound_wave_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to audio_wave (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called
    
    % Hint: place code in OpeningFcn to populate audio_wave

% --- Executes when entered data in editable cell(s) in music_table.
function music_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to music_table (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in table_file_path.
function table_file_path_Callback(hObject, eventdata, handles)
% hObject    handle to table_file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns table_file_path contents as cell array
%        contents{get(hObject,'Value')} returns selected item from table_file_path


% --- Executes during object creation, after setting all properties.
function table_file_path_CreateFcn(hObject, eventdata, handles)
% hObject    handle to table_file_path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function music_table_CreateFcn(hObject, eventdata, handles)
% hObject    handle to music_table (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function tune_Callback(hObject, eventdata, handles)
% hObject    handle to tune (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tune as text
%        str2double(get(hObject,'String')) returns contents of tune as a double


% --- Executes during object creation, after setting all properties.
function tune_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tune (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function length_of_beat_Callback(hObject, eventdata, handles)
% hObject    handle to length_of_beat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of length_of_beat as text
%        str2double(get(hObject,'String')) returns contents of length_of_beat as a double


% --- Executes during object creation, after setting all properties.
function length_of_beat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to length_of_beat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double


% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function title_of_music_Callback(hObject, eventdata, handles)
% hObject    handle to title_of_music (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of title_of_music as text
%        str2double(get(hObject,'String')) returns contents of title_of_music as a double


% --- Executes during object creation, after setting all properties.
function title_of_music_CreateFcn(hObject, eventdata, handles)
% hObject    handle to title_of_music (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function frequency_sampling_Callback(hObject, eventdata, handles)
    % hObject    handle to frequency_sampling (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Hints: get(hObject,'String') returns contents of frequency_sampling as text
    %        str2double(get(hObject,'String')) returns contents of frequency_sampling as a double
    
    
% --- Executes during object creation, after setting all properties.
function frequency_sampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frequency_sampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [song_pitch,song_length] = my_process_music_data_table(file_path)
% 第一列分辨是bass(1)还是alto(2)还是treble(3)还是pause(0)
% 3	7	0.5
% 3	3	0.5
% 3	2	0.75
% 3	3	0.25
% 3	1	1
% 0	0	0.5
% 3	1	0.5
% 2	7	0.5
% 3	1	0.5
% 3	3	0.75
% 3	2	0.25
% 2	6	1
% 0	0	0.5
% 2	3	0.5
% 3	1	0.5
% 2	7	0.75
% 2	6	2
    data=xlsread(file_path);%读取文件datafile.xlsx，并存入data中
    song_pitch = [];
    song_length = data(:,3); % 取出第三列（音长）
    song_pitch_1 = data(:,1);% 取出第一列
    song_pitch_2 = data(:,2);% 取出第一列
    for i = 1:length(song_pitch_1)
        if song_pitch_1(i) == 1
            song_pitch(i) = song_pitch_2 (i)
        elseif song_pitch_1(i) == 2
            song_pitch(i) = song_pitch_2 (i) + 7
        elseif song_pitch_1(i) == 3
            song_pitch(i) = song_pitch_2 (i) + 14
        elseif song_pitch_1(i) == 0
            song_pitch(i) = 22
        else
            error('    Input error!');
        end
    end
    song_pitch = song_pitch';


function tone_sampling = my_play_music(song_pitch,song_length, main_tune,length_of_beat,frequency_sampling,harmonic)

    tunes = get_tunes(main_tune);

    tone_sampling = [];
    [~, padding] = envelope(0);
    last_padding = 0;
    for i = 1 : 1 : length(song_length)
        f = tunes(song_pitch(i));                                                               %The pitch of each tone
        length_of_each_tone = song_length(i) * length_of_beat;                                  %The length of each tone
        length_of_each_tone_padding = frequency_sampling * length_of_each_tone * padding;
        t = linspace(0, length_of_each_tone * padding - 1 / frequency_sampling, length_of_each_tone_padding)';
        tone_sampling_temp = zeros(size(t));
        for j = 1 : 1 : length(harmonic)
            tone_sampling_temp =tone_sampling_temp + harmonic(j) * sin(2 * pi * j * f * t);
        end
        tone_sampling_temp = envelope(t/length_of_each_tone) .* tone_sampling_temp;
        if (last_padding == 0)
            tone_sampling = [tone_sampling; tone_sampling_temp];
        else
            tone_sampling = [tone_sampling(1:end-last_padding); tone_sampling(end-last_padding)+tone_sampling_temp(1:last_padding); tone_sampling_temp(last_padding+1:end)];
        end
        last_padding = round(length_of_each_tone_padding - frequency_sampling * length_of_each_tone);
    end

    % plot([0 : length(tone_sampling) - 1] / frequency_sampling, tone_sampling);
    % title("sound waveform");

    % sound(tone_sampling, frequency_sampling);                       %play sampled tones
    % audiowrite('music_05.wav', tone_sampling, frequency_sampling);  %store .wav
    
function [y, beta_param] = envelope(x)
    alpha_param = 0.075;
    gamma_param = 0.75;
    beta_param = 1.1;
    delta_param = 1 + (gamma_param-alpha_param)./(alpha_param-beta_param);
    A = -1./(alpha_param.^2);
    B = (delta_param-1)./(gamma_param-alpha_param).^2;
    C = (gamma_param-alpha_param)./(gamma_param-beta_param) .* B;

    y = ...
        (x >= 0 & x < alpha_param) .* (A .* (x-alpha_param).^2 + 1) + ...
        (x >= alpha_param & x < gamma_param) .* (B .* (x - alpha_param).^2 + 1) + ...
        (x >= gamma_param & x < beta_param) .* C .* (x - beta_param).^2;
    
function tunes = get_tunes(major)
    tunes = zeros([7 4]);
    ratio = 2^(1/12);

    switch (major)
        case 'A'
            base_A = 1;
        case 'B'
            base_A = 7;
        case 'C'
            base_A = 6;
        case 'D'
            base_A = 5;
        case 'E'
            base_A = 4;
        case 'F'
            base_A = 3;
        case 'G'
            base_A = 2;
        otherwise
            error('Error: The parameter major is invalid!');
    end

    tunes(base_A, 1) = 220;
    tunes(base_A, 2) = 440;
    tunes(base_A, 3) = 880;

    switch base_A
        case 2
            scale_diffs = [-2, 0, 2, 3, 5, 7, 8]';
        case 3
            scale_diffs = [-4, -2, 0, 2, 3, 5, 7]';
        case 4
            scale_diffs = [-5, -4, -2, 0, 2, 3, 5]';
        case 5
            scale_diffs = [-7, -5, -4, -2, 0, 2, 3]';
        case 6
            scale_diffs = [-9, -7, -5, -4, -2, 0, 2]';
        case 7
            scale_diffs = [-10, -9, -7, -5, -4, -2, 0]';
        case 1
            base_A = 1;
            tunes = tunes / 2;
            scale_diffs = [0, 2, 3, 5, 7, 8, 10]';
        otherwise
            error('Error: Unknown error!');
    end

    real_diffs = ratio.^scale_diffs;

    for i = 1 : 7
        tunes(i, 1:end) = tunes(base_A, 1:end) .* real_diffs(i);
    end

function [base_freq_record, component_record] = get_component(wave_path)

    [x, Fs] = audioread(wave_path);
    
    x = mean(x, 2);
    max_Fs = 12000;
    if Fs > max_Fs
        x = resample(x, max_Fs, Fs);
        Fs = max_Fs;
    end
    clear max_Fs;

    % Split music

    disp('Begin to split music...');
    disp('    Processing music...');
    h_wait_bar = waitbar(0, 'Processing music...');

    y1 = abs(x);
    y2WndLen = round(Fs / 10);
    y2 = conv(y1, hanning(y2WndLen));
    y2 = y2(round(y2WndLen/2):end);
    y3 = diff(y2);
    y4 = max(y3, 0);
    y5WndLen = round(Fs / 8);
    y5 = conv(y4, hanning(y5WndLen));
    y5 = y5(round(y5WndLen/2):end);

    figure(1);
    subplot(6, 1, 1);
    plot([0:length(x)-1] / Fs, x);
    title('\itx');
    subplot(6, 1, 2);
    plot([0:length(y1)-1] / Fs, y1);
    title('{\ity}_1');
    subplot(6, 1, 3);
    plot([0:length(y2)-1] / Fs, y2);
    title('{\ity}_2');
    subplot(6, 1, 4);
    plot([0:length(y3)-1] / Fs, y3);
    title('{\ity}_3');
    subplot(6, 1, 5);
    plot([0:length(y4)-1] / Fs, y4);
    title('{\ity}_4');
    subplot(6, 1, 6);
    plot([0:length(y5)-1] / Fs, y5);
    title('{\ity}_5');

    clear y1 y2 y3 y4 y2WndLen y5WndLen;
    waitbar(1, h_wait_bar);
    disp('    OK.');
    disp('    Finding split point...');
    waitbar(0, h_wait_bar, 'Finding split point...');

    dy5 = diff(y5);
    is_positive = (dy5(1) > 0);
    primary_find_idx = [];
    for i = 2 : 1 : length(dy5)
        if ((dy5(i) > 0) ~= is_positive)
            is_positive = ~is_positive;
            if is_positive == false
                primary_find_idx = [primary_find_idx; i];
            end
        end
    end
    waitbar(1 / 3, h_wait_bar);

    min_find_len = 3;
    if length(primary_find_idx) < min_find_len * 2
        error('    Error: The music is too short!');
    end

    primary_find_val = y5(primary_find_idx);
    primary_sort_res = sort(primary_find_val, 'descend');
    level = mean(primary_sort_res(1:min_find_len)) ./ 20;
    secondary_find_idx = primary_find_idx(y5(primary_find_idx) >= level);
    
    waitbar(2 / 3, h_wait_bar);

    thirdary_find_idx = secondary_find_idx;
    min_n_time_inteval = 0.05 * Fs;
    i = 2;
    del_cnt = 0;
    while i <= length(thirdary_find_idx) - del_cnt
        if thirdary_find_idx(i) - thirdary_find_idx(i - 1) < min_n_time_inteval
            if y5(thirdary_find_idx(i - 1)) < y5(thirdary_find_idx(i))
                thirdary_find_idx(i - 1) = thirdary_find_idx(i);
            end
            thirdary_find_idx(i:end-1) = thirdary_find_idx(i+1:end);
            del_cnt = del_cnt + 1;
        else
            i = i + 1;
        end
    end
    thirdary_find_idx = thirdary_find_idx(1:end-del_cnt);

    find_idx = thirdary_find_idx;

    subplot(6, 1, 6);
    hold on
    plot((find_idx-1)/Fs, y5(find_idx), 'ro');
    subplot(6, 1, 1);
    hold on
    plot((find_idx-1)/Fs, zeros([length(find_idx), 1]), 'yo');

    clear y5 dy5 primary_find_val primary_sort_res primary_find_idx secondary_find_idx;
    waitbar(1, h_wait_bar);
    close(h_wait_bar);
    disp('    OK.');
    disp('Spliting music finished!');
    disp('Begin to find basic frequency...');
    disp('    Preparing...');
    h_wait_bar = waitbar(0, 'Preparing for finding...');

    amp_tolerant_rate_of_base = 0.1;

    base_freqs = [];
    base_freq_idxs = [];
    Xs = cell([0 0]);
    T1s = [];

    disp('    OK.');
    disp('    Finding...');
    waitbar(0, h_wait_bar, 'Finding basic frequencies...');

    % Find basic freqency in each piece of music

    for i = 2 : 1 : length(find_idx)

        this_b = find_idx(i - 1);
        this_e = find_idx(i);
        this_x = x(this_b : this_e);
        this_repeat = 10;
        for j = 1 : 1 : this_repeat
            this_x = [this_x; this_x];
        end
        this_repeat = 2^this_repeat;
        this_X = abs(fft(this_x));
        T1 = (this_e - this_b + 1) * this_repeat / Fs;

        dc_comp = this_X(1);  % Exclude DC component
        this_X = this_X(2:end);
        [max_val, max_idx] = max(this_X);
        amp_level = max_val / 3;
        over_level = (this_X > amp_level);
        [~, first_idx] = max(over_level);

        if first_idx == max_idx
            base_freq_idx = max_idx;
        else

            candidate_idx = 0;
            candidate_times = 0;
            candidate_times_ret = 0;

            while first_idx < max_idx

                % disp("* " + first_idx/T1);

                tmp_times = max_idx / first_idx;

                % disp(tmp_times);

                tmp_ret = abs(tmp_times - round(tmp_times));
                if tmp_times < 20 && tmp_ret < 0.1
                    if (tmp_times < candidate_times && candidate_idx ~= 0)
                        break;
                    end
                    if (candidate_idx == 0 || tmp_ret < candidate_times_ret)
                        should_be_candidate = true;
                        %{
                        if tmp_times ~= 1
                            time_of_this_over_times = 0;
                            time_of_max_over_times = 0;
                            is_times = @(p, q) abs((p/q) - round(p/q)) < 0.1;
                            for k = 1 : 1 : 15 % amp_tolerant_rate_of_base
                                time_of_this = k * first_idx;
                                tolerant_idxs = round(first_idx * amp_tolerant_rate_of_base);
                                if is_times(time_of_this, max_idx)
                                    if max(this_X(time_of_this-tolerant_idxs : time_of_this+tolerant_idxs)) > amp_level
                                        time_of_max_over_times = time_of_max_over_times + 1;
                                    end
                                else
                                    if max(this_X(time_of_this-tolerant_idxs : time_of_this+tolerant_idxs)) > amp_level
                                        time_of_this_over_times = time_of_this_over_times + 1;
                                    end
                                end
                            end

                            if time_of_max_over_times >= round(time_of_this_over_times*time_of_this_over_times*0.5-time_of_this_over_times*0.5+3)-0.05;
                                hould_be_candidate = false;
                            else 
                            end
                        end
                        %}
                        if should_be_candidate == true
                            candidate_idx = first_idx;
                            candidate_times = round(tmp_times);
                            candidate_times_ret = tmp_ret;
                        end
                    end
                end
                [val, tmp_first_idx] = max(over_level(first_idx + 1 : end));
                first_idx = tmp_first_idx + first_idx;
            end

            if candidate_idx >= max_idx || candidate_times == 1 || candidate_idx == 0
                base_freq_idx = max_idx;
            else
                base_freq_idx = candidate_idx;
            end
        end

        % disp("    Origin: " + (i-1) + " ~ " + i + ": " + (base_freq_idx / T1));

        base_freqs = [base_freqs; base_freq_idx / T1];
        base_freq_idxs = [base_freq_idxs; base_freq_idx];
        Xs = [Xs; this_X];
        T1s = [T1s; T1];

        % sound(this_x(1:min(Fs, end)), Fs);

        %{
        figure(2);
        this_X = [dc_comp; this_X];
        plot([0:length(this_X)-1] * (1/T1), this_X);
        disp((i - 1) + " ~ " + i + ": ");
        pause
        %}
        
        waitbar(i / length(find_idx), h_wait_bar);
    end

    clear this_b this_e this_x this_X this_repeat;
    clear over_level candidate_idx candidate_times candidate_times_ret;
    waitbar(1, h_wait_bar);
    disp('    OK.');
    disp('    Checking legitimacy of frequencies found...');
    waitbar(0, h_wait_bar, 'Checking legitimacy...');

    min_freq = 110;
    max_freq = 1500;
    freq_too_low = base_freqs < min_freq;
    freq_too_high = base_freqs > max_freq;
    legal_freq = ~(freq_too_low | freq_too_high);
    if sum(freq_too_low | freq_too_high) ~= 0
        disp('    Warning: A basic frequency is not in the list and will be discarded!');
    end
    base_freqs = base_freqs(legal_freq);
    base_freq_idxs = base_freq_idxs(legal_freq);
    Xs = Xs(legal_freq);
    T1s = T1s(legal_freq);
    std_freqs = generate_std_freqs(min_freq, max_freq);
    modified_base_freqs = nearest_search(std_freqs, base_freqs);

    clear freq_too_low freq_too_high legal_greq;
    waitbar(1, h_wait_bar);
    disp('    OK.');
    disp('Finding basic frequency finished!');
    disp('Begin to get components...');
    disp('    Preparing...');
    waitbar(0, h_wait_bar, 'Preparing for getting components...');

    % Sort by freqs

    [modified_base_freqs, sort_freq_idx] = sort(modified_base_freqs);
    base_freqs = base_freqs(sort_freq_idx);
    base_freq_idxs = base_freq_idxs(sort_freq_idx);
    Xs = Xs(sort_freq_idx);
    T1s = T1s(sort_freq_idx);

    % Get parameters

    component_record = cell(length(std_freqs), 2); % Column 1: freq; column 2: components
    component_record(:, 1) = num2cell(std_freqs);
    component_record_itr = 1;

    is_float_equal = @(f, g) abs(f - g) / g < 1e-5;

    clear sort_freq_idx;
    disp('    OK.');
    disp('    Getting...');
    waitbar(0, h_wait_bar, 'Getting components...');

    for i = 1 : 1 : length(modified_base_freqs)
        max_times = floor(length(Xs{i}) / 2 / base_freq_idxs(i));
        k = [1 : 1 : max_times]';
        mid_freqs = base_freq_idxs(i) * k;
        tolerant_idxs = round(base_freq_idxs(i) * amp_tolerant_rate_of_base);
        left_freqs = mid_freqs - tolerant_idxs;
        right_freqs = mid_freqs + tolerant_idxs;

        component_res = zeros([max_times, 1]);
        for j = 1 : 1 : max_times
            component_res(j) = max(Xs{i}(left_freqs(j):right_freqs(j)));
        end
        component_res  = component_res / component_res(1);

        while is_float_equal(modified_base_freqs(i), component_record{component_record_itr, 1}) == false
            component_record_itr = component_record_itr + 1;
            if component_record_itr > length(std_freqs)
                error('    Error: unknown error!');
            end
        end

        tmp_org = component_record{component_record_itr, 2};
        % disp('=============');
        % disp(tmp_org);
        if size(tmp_org, 2) == 0
            component_record{component_record_itr, 2} = component_res;
        else
            if size(tmp_org, 1) < length(component_res)
                tmp_org = [tmp_org; zeros([ length(component_res) - size(tmp_org, 1), size(tmp_org, 2)])];
            elseif size(tmp_org, 1) > length(component_res)
                component_res = [component_res; zeros([size(tmp_org, 1) - length(component_res), 1])];
            end
            component_record{component_record_itr, 2} = [tmp_org, component_res];
        end
        % disp('-----');
        % disp(component_record{component_record_itr, 2});
        waitbar(i / length(modified_base_freqs), h_wait_bar);
    end

    clear Xs;
    waitbar(1, h_wait_bar);
    disp('    OK.');
    disp('    Cleaning...');
    waitbar(0, h_wait_bar, 'Cleaning...');

    legal_component = logical(zeros([length(std_freqs), 1]));
    for i = 1 : 1 : length(std_freqs)
        if size(component_record{i, 2}, 2) ~= 0
            legal_component(i) = true;
        end
        waitbar(i / length(std_freqs), h_wait_bar);
    end
    component_record = component_record(legal_component, :);

    %{
    for i = 1 : 1 : 11
        figure(i + 2);
        hold on
        tmp = component_record{i,2};
        for j = 1 : 1 : size(tmp, 2)
            plot(component_record{i,2}(:, j));
        end
        title(string(component_record{i,1}) + " Hz");
    end
    %}

    clear legal_component;
    waitbar(1, h_wait_bar);
    disp('    OK.');
    disp('    Merging multidata...');
    waitbar(0, h_wait_bar, 'Merging multidata...');

    base_freq_record = cell2mat(component_record(:, 1));
    component_record = component_record(:, 2);

    for i = 1 : 1 : length(component_record)
        component_record{i} = mean(component_record{i}, 2);
        waitbar(i / length(component_record), h_wait_bar);
    end

    waitbar(1, h_wait_bar);
    waitbar(1, h_wait_bar, 'Finished!');
    disp('    OK.');
    disp('Getting components finished!');
    close(h_wait_bar);

function std_freqs = generate_std_freqs(base_freq, max_freq)

    std_freqs = [];
    freq_itr = base_freq;
    ratio = 2^(1/12);
    while freq_itr <= max_freq
        std_freqs = [std_freqs; freq_itr];
        freq_itr = freq_itr * ratio;
    end

function [val, idx] = nearest_search(list, target, is_sorted_list)
% nearest_search: find the nearest element

    if nargin < 3
        is_sorted_list = true;
    elseif is_sorted_list == false
        if issorted(list) == false
            list = sort(list);
        end
    end

    if length(list) == 0
       error('The list is empty!');
    end

    if length(target) == 0
        val = [];
        idx = [];
    else
        idx = zeros(size(target));
        for i = 1 : 1 : length(target)
            [~, idx(i)] = min(abs(list - target(i)));
        end
        val = list(idx);
    end
