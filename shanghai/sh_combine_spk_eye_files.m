function flog = sh_combine_spk_eye_files(spk_file,log_file,eye_file,eye_refresh,eye_format,eye_threshold,if_plot)
% function flog = sh_combine_spk_eye_files(spk_file,log_file,eye_file,eye_refresh,eye_format,eye_threshold,if_plot)
% sh_combine_spk_eye_files - combine f21 spike and eye data
%
% INPUTS:  spk_file         - f21 spike file
%          log_file         - f21 log file
%          eye_file         - f21 eye file 
%          eye_refresh      - if re-analyze eye data (def.0), 0 or 1
%          eye_format       - eye data format (def. 0), 0(text) or 1(binary)
%          eye_threshold    - eye movement threshold in Kurtosis (def. 3, normal distribution) 
%          if_plot          - if plot (def.0), 0 (no) or 1 (yes)
%
% OUTPUTS: flog - eye info
%
% Yuxi on 12/21/2004
if nargin<7    if_plot = 0;         end
if nargin<6    eye_threshold = 3;   end % normal distribution
if nargin<5    eye_format = 0;      end
if nargin<4    eye_refresh = 0;     end
if nargin<3
    tmp = length(spk_file);
    eye_file = [spk_file(1:tmp-3) 'eye'];
end
if nargin<2
    tmp = length(spk_file);
    log_file = [spk_file(1:tmp-3) 'log'];
end
flog.spk_file = spk_file;
flog.log_file = log_file;
flog.eye_file = eye_file;
flog.eye_refresh = eye_refresh;
flog.eye_format = eye_format;
flog.eye_threshold = eye_threshold;
flog.if_plot = if_plot;

% retrieve refresh rate
log_lines = fload_log(log_file);
str_temp = fsearch_log(log_lines,'TestInfo','RefreshRate','1');
refresh_rate = sscanf(str_temp,'%f');
unit_per_refresh = 10000/refresh_rate;
flog.refresh_rate = refresh_rate;
flog.unit_per_refresh = unit_per_refresh;

% load eye data
eye_MAT_file = [eye_file '.MAT'];
if eye_refresh~=0 | ~exist(eye_MAT_file,'file') 
    eye_flog = sh_ana_eye_file(eye_file,eye_format,eye_threshold,if_plot);
else
    tmp = load(eye_MAT_file);
    eye_flog = tmp.flog;
end
eye_time_ms = eye_flog.eye_data.eye_data(:,1);
eye_time_num = length(eye_time_ms);
eye_data_good = eye_flog.eye_data_good;
flog.eye_flog = eye_flog;
flog.eye_time_ms = eye_time_ms;
flog.eye_time_num = eye_time_num;
flog.eye_data_good = eye_data_good;

% spike file & histogram
spikes = load_spk_file(spk_file);
max_ms = max(max(eye_time_ms),max(spikes)*0.1);
num_bins = ceil(max_ms/(unit_per_refresh*0.1));
tmp = y_bin_3(spikes,unit_per_refresh,num_bins);
histogram = tmp.histogram;
flog.spikes = spikes;
flog.num_bins = num_bins;
flog.histogram = histogram;

% eye file & histogram
dms = zeros(eye_time_num-1,1);
for i=1:eye_time_num-1
    dms(i) = eye_time_ms(i+1)-eye_time_ms(i);
end
if (max(dms(:))-min(dms(:)))>1 % they should be identical
    fprintf('ERROR: measurement interval in eye file %s is not consistent.\n',eye_file);
    return;
end
dms_unit = round(max(dms)*10);
eye_virtual_spike_array = zeros(eye_time_num*dms_unit,1);
for i=1:eye_time_num
    p1 = round(max(eye_time_ms(i)*10-dms_unit*0.5,1));
    p2 = round(eye_time_ms(i)*10+dms_unit*0.5);
    eye_virtual_spike_array(p1:p2) = eye_data_good(i);
end
eye_virtual_spikes = find(eye_virtual_spike_array==1);
eye_tmp = y_bin_3(eye_virtual_spikes,unit_per_refresh,num_bins);
eye_histogram = eye_tmp.histogram;
eye_fixed = eye_histogram>(unit_per_refresh-1);
flog.dms_unit = dms_unit;
%flog.eye_virtual_spikes = eye_virtual_spikes;
flog.eye_histogram = eye_histogram;
flog.eye_fixed = eye_fixed;

% eye stability distribution - in
in_idx = find(eye_fixed==1);
in_tmp = sh_count_continuity(in_idx);
in_counter = in_tmp.counter;
in_len_data = in_tmp.len_data;
flog.in_counter = in_counter;
flog.in_len_data = in_len_data;

% eye stability distribution - out
out_idx = find(eye_fixed==0);
out_tmp = sh_count_continuity(out_idx);
out_counter = out_tmp.counter;
out_len_data = out_tmp.len_data;
flog.out_counter = out_counter;
flog.out_len_data = out_len_data;

% combined in & out
in_out_counter = in_counter + out_counter;
in_out_len_data = zeros(in_out_counter,1);
if eye_fixed(1)==1
    in_out_len_data(1:2:in_out_counter) = in_len_data;
    in_out_len_data(2:2:in_out_counter) = -out_len_data;
else
    in_out_len_data(2:2:in_out_counter) = in_len_data;
    in_out_len_data(1:2:in_out_counter) = -out_len_data;
end
flog.in_out_counter = in_out_counter;
flog.in_out_len_data = in_out_len_data;

% save MAT file
save([spk_file '.MAT'],'flog');

% figure
if if_plot==0    return;   end

%% eye position versus time
fig=figure('NumberTitle','off','Name',spk_file,'Position',[100 250 600 400]);
hold on

max_tmp = max(abs(in_out_len_data));
min_tmp = min(abs(in_out_len_data));
l = 0;
r = in_out_counter;
b = min_tmp - (max_tmp - min_tmp)*0.15;
t = max_tmp + (max_tmp - min_tmp)*0.15;

set(gca,'Position',[0.12 0.10 0.55 0.80],'Box','On');
xlabel('Sessions');
ylabel('<- Out of focus (# of frames) In focus ->');
tmp = spk_file;
tmp(tmp=='_') = '-';
title(tmp);
axis([l r b t]);

%plot(1:in_out_counter,in_out_len_data);
for i=1:in_out_counter
    if in_out_len_data(i)>0
        plot([i i],[0 in_out_len_data(i)],'b','LineWidth',2);
    else
        plot([i i],[0 in_out_len_data(i)],'r','LineWidth',2);
    end
end

grid on

% text info 
txt_info{1} = sprintf('Refresh Rate = %0.4f Hz',refresh_rate);
txt_info{2} = sprintf('Frame Size = %0.3f ms',1000/refresh_rate);
text(r+(r-l)*0.05,t-(t-b)*0.0,txt_info,'VerticalAlignment','top','HorizontalAlignment','left');

% exit
return;
