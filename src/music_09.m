clear;
close all;
clc;
% Plot
% function my_plot(x,y1,y2,y3,y4,y5,frequency_sampling,find_idx,component_record)
% load_file
% function [x,frequency_sampling] = my_load_file(file_path)
% Split music
% function [y1,y2,y3,y4,y5] = my_split_music(x,frequency_sampling)
% Find idx
% function [find_idx] = my_find_idx(y5,frequency_sampling)
% Find basic freqency in each piece of music
% function [base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs,std_freqs] = my_find_base_freq(find_idx,frequency_sampling,x)
% Sort by freqs
% function [base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs] = my_sort_by_freqs(base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs)
% Get parameters
% function [component_res,component_record] = my_get_parameters(std_freqs,modified_base_freqs,base_freq_idxs,Xs)
% record
% function component_record= my_record(std_freqs,component_record)
% generate_std_freqs
% function std_freqs = generate_std_freqs(base_freq, max_freq)
% nearest_search: find the nearest element
% function [val, idx] = nearest_search(list, target, is_sorted_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,frequency_sampling] = my_load_file('resource/fmt.wav');
[y1,y2,y3,y4,y5] =  my_split_music(x,frequency_sampling);
[find_idx] = my_find_idx(y5,frequency_sampling);
[base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs,std_freqs] = my_find_base_freq(find_idx,frequency_sampling,x);
[base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs] = my_sort_by_freqs(base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs);
[component_res,component_record] = my_get_parameters(std_freqs,modified_base_freqs,base_freq_idxs,Xs);
[component_record,base_freq_record] = my_record(std_freqs,component_record);
my_plot(x,y1,y2,y3,y4,y5,frequency_sampling,find_idx,component_record);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot
function my_plot(x,y1,y2,y3,y4,y5,frequency_sampling,find_idx,component_record)
    % Automatically time the beginning and end of notes
    figure(1);
    subplot(6, 1, 1);
    plot([0:length(x)-1] / frequency_sampling, x);
    title('x-原信号');
    subplot(6, 1, 2);
    plot([0:length(y1)-1] / frequency_sampling, y1);
    title('{y}_1-位移绝对值');
    subplot(6, 1, 3);
    plot([0:length(y2)-1] / frequency_sampling, y2);
    title('{y}_2-能量');
    subplot(6, 1, 4);
    plot([0:length(y3)-1] / frequency_sampling, y3);
    title('{y}_3-变化率');
    subplot(6, 1, 5);
    plot([0:length(y4)-1] / frequency_sampling, y4);
    title('{y}_4-半波整流');
    subplot(6, 1, 6);
    plot([0:length(y5)-1] / frequency_sampling, y5);
    title('{y}_5-加窗平滑');

    %Extraction of harmonic components
    subplot(6, 1, 6);
    hold on
    plot((find_idx-1)/frequency_sampling, y5(find_idx), 'go');
    subplot(6, 1, 1);
    hold on
    plot((find_idx-1)/frequency_sampling, zeros([length(find_idx), 1]), 'ro');

    figure(2);
    for i = 1 : 1 : 11
        subplot(4,3,i);
        hold on
        tmp = component_record{i,2};
        for j = 1 : 1 : size(tmp, 2)
            plot(component_record{i,2}(:, j));
        end
        title(string(component_record{i,1}) + " Hz");
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load_file
function [x,frequency_sampling] = my_load_file(file_path)
    [x, frequency_sampling] = audioread(file_path);
    x = mean(x, 2);
    max_frequency_sampling = 12000;
    if frequency_sampling > max_frequency_sampling
        x = resample(x, max_frequency_sampling, frequency_sampling);
        frequency_sampling = max_frequency_sampling;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split music
function [y1,y2,y3,y4,y5] = my_split_music(x,frequency_sampling)
    y1 = abs(x);
    y2WndLen = round(frequency_sampling / 10);
    y2 = conv(y1, hanning(y2WndLen));
    y2 = y2(round(y2WndLen/2):end);
    y3 = diff(y2);
    y4 = max(y3, 0);
    y5WndLen = round(frequency_sampling / 8);
    y5 = conv(y4, hanning(y5WndLen));
    y5 = y5(round(y5WndLen/2):end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find idx
function [find_idx] = my_find_idx(y5,frequency_sampling)
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
    min_find_len = 3;
    if length(primary_find_idx) < min_find_len * 2
        error('    Error: The music is too short!');
    end
    primary_find_val = y5(primary_find_idx);
    primary_sort_res = sort(primary_find_val, 'descend');
    level = mean(primary_sort_res(1:min_find_len)) ./ 20;
    secondary_find_idx = primary_find_idx(y5(primary_find_idx) >= level);
    thirdary_find_idx = secondary_find_idx;
    min_n_time_inteval = 0.05 * frequency_sampling;
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find basic freqency in each piece of music
function [base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs,std_freqs] = my_find_base_freq(find_idx,frequency_sampling,x)
    base_freqs = [];
    base_freq_idxs = [];
    Xs = cell([0 0]);
    T1s = [];
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
        T1 = (this_e - this_b + 1) * this_repeat / frequency_sampling;
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
                tmp_times = max_idx / first_idx;
                tmp_ret = abs(tmp_times - round(tmp_times));
                if tmp_times < 20 && tmp_ret < 0.1
                    if (tmp_times < candidate_times && candidate_idx ~= 0)
                        break;
                    end
                    if (candidate_idx == 0 || tmp_ret < candidate_times_ret)
                        should_be_candidate = true;
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
        base_freqs = [base_freqs; base_freq_idx / T1];
        base_freq_idxs = [base_freq_idxs; base_freq_idx];
        Xs = [Xs; this_X];
        T1s = [T1s; T1];
    end
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort by freqs
function [base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs] = my_sort_by_freqs(base_freqs,base_freq_idxs,Xs,T1s,modified_base_freqs)
    [modified_base_freqs, sort_freq_idx] = sort(modified_base_freqs);
    base_freqs = base_freqs(sort_freq_idx);
    base_freq_idxs = base_freq_idxs(sort_freq_idx);
    Xs = Xs(sort_freq_idx);
    T1s = T1s(sort_freq_idx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get parameters
function [component_res,component_record] = my_get_parameters(std_freqs,modified_base_freqs,base_freq_idxs,Xs)
    component_record = cell(length(std_freqs), 2); % Column 1: freq; column 2: components
    component_record(:, 1) = num2cell(std_freqs);
    component_record_itr = 1;
    is_float_equal = @(f, g) abs(f - g) / g < 1e-5;
    amp_tolerant_rate_of_base = 0.1;
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
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% record
function [component_record,base_freq_record]= my_record(std_freqs,component_record)
    legal_component = logical(zeros([length(std_freqs), 1]));
    for i = 1 : 1 : length(std_freqs)
        if size(component_record{i, 2}, 2) ~= 0
            legal_component(i) = true;
        end
    end
    component_record = component_record(legal_component, :);

    base_freq_record = cell2mat(component_record(:, 1));
    component_record_temp = component_record(:, 2);
    for i = 1 : 1 : length(component_record_temp)
        component_record_temp{i} = mean(component_record_temp{i}, 2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_std_freqs
function std_freqs = generate_std_freqs(base_freq, max_freq)
    std_freqs = [];
    freq_itr = base_freq;
    ratio = 2^(1/12);
    while freq_itr <= max_freq
        std_freqs = [std_freqs; freq_itr];
        freq_itr = freq_itr * ratio;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nearest_search: find the nearest element
function [val, idx] = nearest_search(list, target, is_sorted_list)
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
end