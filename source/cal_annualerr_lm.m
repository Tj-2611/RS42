function [diff_list] = cal_annualerr_lm(x, optV_tmp, model, method, spx_cell, vix_cell, total_num)
% calculating the annual errors using levenberg-marquardt algorithm.
t_num = size(optV_tmp, 2);
diff_cell = cell(t_num, 1);
if isrow(x)
    x = x';
end
parfor i = 1:t_num
    params = [optV_tmp(:, i); x];
    spx_tmp = spx_cell{i};
    vix_tmp = vix_cell{i};
    diff_cell{i} = cal_volsurfsep(params, spx_tmp, vix_tmp, model, method, 'lm');
end
diff_list = zeros(total_num, 1)*nan;
count = 1;
for i = 1:t_num
    tmp = diff_cell{i};
    m = length(tmp);
    diff_list(count:count+m-1, 1) = tmp/sqrt(total_num);
    count = count + m;
end
end