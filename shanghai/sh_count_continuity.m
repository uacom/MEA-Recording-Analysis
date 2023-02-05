function flog = sh_count_continuity(idx)
% function flog = sh_count_continuity(idx)
%
% INPUTS    data - idx vector, like [1 2 4 6 7 8]
%
% OUTPUTS flog.
%
% Yuxi Fu 12/25/2004

idx_num = length(idx);
p1 = 1;
counter = 0;
len_data = zeros(idx_num,1);
for i=2:idx_num
    if idx(i)-idx(i-1)>1
        p2 = i-1;
        cur_len = p2-p1+1;
        counter = counter + 1;
        len_data(counter) = cur_len;
        p1 = i;
    else
        ;
    end
end
% last
p2 = idx_num;
cur_len = p2-p1+1;
counter = counter + 1;
len_data(counter) = cur_len;
len_data = len_data(1:counter);

% flog
flog.idx = idx;
flog.idx_num = idx_num;
flog.counter = counter;
flog.len_data = len_data;

% exit
return;
