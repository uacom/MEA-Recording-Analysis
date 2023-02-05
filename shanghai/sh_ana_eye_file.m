function flog = sh_ana_eye_file(eye_file,eye_format,eye_threshold,if_plot)
% function flog = sh_ana_eye_file(eye_file,eye_format,eye_threshold,if_plot)
% sh_ana_eye_file - analyze f21 eye data file
%
% INPUTS:  eye_file         - f21 eye file 
%          eye_format       - eye data format (def. 0), 0(text) or 1(binary)
%          eye_threshold    - eye movement threshold in Kurtosis (def. 3, normal distribution) 
%          if_plot          - if plot (def.0), 0 (no) or 1 (yes)
%
% OUTPUTS:  flog - eye info
%           flog.eye_data
%           flog.eye_x_mean
%           flog.eye_x_threshold
%           flog.eye_x_kurt
%           flog.eye_x_percentage
%           flog.eye_y_mean
%           flog.eye_y_threshold
%           flog.eye_y_kurt
%           flog.eye_y_percentage
%           flog.eye_data_good
%           flog.eye_data_percentage
%
% Yuxi on 12/21/2004

if nargin<4    if_plot = 0;         end
if nargin<3    eye_threshold = 3;   end % normal distribution
if nargin<2    eye_format = 0;      end
flog.eye_file = eye_file;
flog.eye_format = eye_format;
flog.eye_threshold = eye_threshold;
flog.if_plot = if_plot;

% constants
EYE_X_IDX = 2; % left eye (right eye, use 4)
EYE_Y_IDX = 3; % left eye (right eye, use 5)
KURTOSIS_STEP_SIZE = 0.8;

% load eye data
if eye_format==0
    eye_data = sh_load_txt_eye_file(eye_file);
else
    eye_data = sh_load_bin_eye_file(eye_file);
end
flog.eye_data = eye_data;

% find X threshold
eye_x_dat = eye_data.eye_data(:,EYE_X_IDX);
eye_x_tmp = sh_find_kurtosis_threshold(eye_x_dat,eye_threshold,KURTOSIS_STEP_SIZE);
eye_x_mean = eye_x_tmp.data_mean;
eye_x_threshold = eye_x_tmp.data_threshold;
eye_x_kurt = eye_x_tmp.data_kurt;
eye_x_percentage = eye_x_tmp.data_percentage;
flog.eye_x_mean = eye_x_mean;
flog.eye_x_threshold = eye_x_threshold;
flog.eye_x_kurt = eye_x_kurt;
flog.eye_x_percentage = eye_x_percentage;

% find Y threshold
eye_y_dat = eye_data.eye_data(:,EYE_Y_IDX);
eye_y_tmp = sh_find_kurtosis_threshold(eye_y_dat,eye_threshold,KURTOSIS_STEP_SIZE);
eye_y_mean = eye_y_tmp.data_mean;
eye_y_threshold = eye_y_tmp.data_threshold;
eye_y_kurt = eye_y_tmp.data_kurt;
eye_y_percentage = eye_y_tmp.data_percentage;
flog.eye_y_mean = eye_y_mean;
flog.eye_y_threshold = eye_y_threshold;
flog.eye_y_kurt = eye_y_kurt;
flog.eye_y_percentage = eye_y_percentage;

% combine X & Y
num_eye_data = eye_data.num_dat_lines;
eye_data_good = ones(num_eye_data,1);
for i=1:num_eye_data
    if eye_x_dat(i)<eye_x_mean-eye_x_threshold | eye_x_dat(i)>eye_x_mean+eye_x_threshold
        eye_data_good(i) = 0;
    end
    if eye_y_dat(i)<eye_y_mean-eye_y_threshold | eye_y_dat(i)>eye_y_mean+eye_y_threshold
        eye_data_good(i) = 0;
    end
end
flog.eye_data_good = eye_data_good;
flog.eye_data_percentage = sum(eye_data_good)/num_eye_data*100;

% save MAT file
save([eye_file '.MAT'],'flog');

% figure
if if_plot==0    return; end

%% eye position versus time
fig=figure('NumberTitle','off','Name',eye_file,'Position',[100 250 500 400]);
hold on

time_sec = eye_data.eye_data(:,1)/1000; % time in second

subplot(2,1,1);
hold on
l = 0;
r = max(time_sec);
tmpmin = min(eye_x_dat);
tmpmax = max(eye_x_dat);
b = tmpmin - (tmpmax-tmpmin)*0.05;
t = tmpmax + (tmpmax-tmpmin)*0.05;
ylabel('Eye X (cm)');
tmp = eye_file;
tmp(tmp=='_') = '-';
title(tmp);
axis([l r b t]);
plot(time_sec,eye_x_dat,'b');
grid on

subplot(2,1,2);
hold on
l = 0;
r = max(time_sec);
tmpmin = min(eye_y_dat);
tmpmax = max(eye_y_dat);
b = tmpmin - (tmpmax-tmpmin)*0.05;
t = tmpmax + (tmpmax-tmpmin)*0.05;
xlabel('Time (sec)');
ylabel('Eye Y (cm)');
axis([l r b t]);
plot(time_sec,eye_y_dat,'b');
grid on

%% eye scatter plot
fig=figure('NumberTitle','off','Name',eye_file,'Position',[300 300 600 400]);
hold on

eye_box_range = 3;
tmp = max(eye_x_threshold,eye_y_threshold);
l = eye_x_mean - tmp*eye_box_range;
r = eye_x_mean + tmp*eye_box_range;
b = eye_y_mean - tmp*eye_box_range;
t = eye_y_mean + tmp*eye_box_range;

set(gca,'Position',[0.08 0.10 0.55 0.80],'Box','On');
xlabel('Eye X (cm)');
ylabel('Eye Y (cm)');
tmp = eye_file;
tmp(tmp=='_') = '-';
title(tmp);
axis([l r b t]);

% in box points
in_idx = find(eye_data_good==1);
in_x = eye_x_dat(in_idx);
in_y = eye_y_dat(in_idx);
plot(in_x,in_y,'*b','LineStyle','none');

% out box points
out_idx = find(eye_data_good==0);
out_x = eye_x_dat(out_idx);
out_y = eye_y_dat(out_idx);
plot(out_x,out_y,'*k','LineStyle','none');

% center
plot([eye_x_mean eye_x_mean],[eye_y_mean eye_y_mean],'or','LineStyle','none');

% text info 
txt_info{1} = sprintf('Kurtosis Thres. = %6.3f',eye_threshold);
txt_info{2} = '';
txt_info{3} = sprintf('X Mean = %6.3f',eye_x_mean);
txt_info{4} = sprintf('X Threshold = %6.3f',eye_x_threshold);
txt_info{5} = sprintf('X Kurtosis = %6.3f',eye_x_kurt);
txt_info{6} = sprintf('X %% = %6.3f',eye_x_percentage);
txt_info{7} = '';
txt_info{8} = sprintf('Y Mean = %6.3f',eye_y_mean);
txt_info{9} = sprintf('Y Threshold = %6.3f',eye_y_threshold);
txt_info{10} = sprintf('Y Kurtosis = %6.3f',eye_y_kurt);
txt_info{11} = sprintf('Y %% = %6.3f',eye_y_percentage);
txt_info{12} = '';
txt_info{13} = sprintf('X&Y %% = %6.3f',flog.eye_data_percentage);
text(r+(r-l)*0.1,t-(t-b)*0.0,txt_info,'VerticalAlignment','top','HorizontalAlignment','left');

grid on

% exit
return;
