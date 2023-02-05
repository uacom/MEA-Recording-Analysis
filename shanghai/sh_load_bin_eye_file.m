function flog = sh_load_bin_eye_file(eye_file)
% function flog = sh_load_bin_eye_file(eye_file)
% sh_load_txt_eye_file - LOADING F21 EYE DATA FILE IN BINARY FORMAT
%
% INPUTS:  eye_file - eye file in BINARY format
%
% OUTPUTS: flog - eye info
%          flog.file_id: 
%          flog.file_name:
%          flog.num_dat_lines: 
%          flog.eye_data: (num_dat_lines,5) - TIME(ms), LH, LV, RH, RV
%          flog.IsFailed: 
%
% Yuxi on 12/20/2004

flog.file_id = ''; 
flog.file_name = '';
flog.num_dat_lines = 0; 
flog.eye_data = [];

% exit
flog.IsFailed = 1;
return;
