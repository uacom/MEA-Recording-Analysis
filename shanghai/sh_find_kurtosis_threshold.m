function flog = sh_find_kurtosis_threshold(data,kurt_threshold,find_step)
%function flog = sh_find_kurtosis_threshold(data,kurt_threshold,find_step)
%
% INPUTS    data            - data
%           kurt_threshold  - kurtosis threshold (def.3, normal distribution)
%           find_step       - step of threshold change (def. 0.9), (0,1)
%   
% OUTPUTS   flog.data_final: 
%           flog.data_mean: 
%           flog.data_threshold:
%           flog.data_kurt: 
%           flog.data_percentage:
%
% Yuxi Fu, 12/22/2004

if nargin<3    find_step = 0.9;     end
if nargin<2    kurt_threshold = 3;  end % normal distribution
flog.data_ori = data;
flog.kurt_threshold = kurt_threshold;
flog.find_step = find_step;

% search
data_num = length(data);
data_threshold = max(abs(data-mean(data)));
data_kurt = 1000;
while data_kurt>kurt_threshold
    data_threshold = data_threshold*find_step;
    low = mean(data) - data_threshold;
    high = mean(data) + data_threshold;
    data = data(data>=low & data<=high);
    data_kurt = kurtosis(data);
    fprintf('Kurtosis Search: mean=%6.2f threshold=%6.2f kurt=%6.2f percentage=%6.2f\n',mean(data),...
        data_threshold,data_kurt,length(data)/data_num*100);
end
flog.data_final = data;
flog.data_mean= mean(data);
flog.data_threshold = data_threshold;
flog.data_kurt = data_kurt;
flog.data_percentage = length(data)/data_num*100;

% exit
return;