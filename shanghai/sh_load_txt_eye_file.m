function flog = sh_load_txt_eye_file(eye_file)
% function flog = sh_load_txt_eye_file(eye_file)
% sh_load_txt_eye_file - LOADING F21 EYE DATA FILE IN TEXT FORMAT
%
% INPUTS:  eye_file - eye file in TXT format
%
% OUTPUTS: flog - eye info
%          flog.file_id: 
%          flog.file_name:
%          flog.num_dat_lines: 
%          flog.eye_data: (num_dat_lines,5) - TIME(ms), LH, LV, RH, RV
%          flog.IsFailed: 
%
% Yuxi on 12/20/2004

fprintf('Loading eye data file ''%s'' ... ',eye_file);
lines = fload_log(eye_file);
fprintf('finished\n');

file_id = lines{1};
file_name = lines{2};
flog.file_id = file_id;
flog.file_name = file_name;

if ~strcmp(file_id,'f21mlv.vm11.vi: Eye Tracking Data')
    fprintf('ERROR: file ''%s'' is not supported.\n',eye_file);
    flog.IsFailed = 1;
    return;
end

% number of data lines
num_lines_p0 = 23;
str_tmp = lines{3};
num_dat_lines = str2num(str_tmp(num_lines_p0:length(str_tmp)));
flog.num_dat_lines = num_dat_lines;

% skip one lines - captions
%   Time(ms)        LX        LY        RX        RY

% data
eye_data = zeros(num_dat_lines,5);
for i=1:num_dat_lines
    eye_data(i,:) = str2num(lines{4+i});
end
flog.eye_data = eye_data;

% exit
flog.IsFailed = 0;
return;
