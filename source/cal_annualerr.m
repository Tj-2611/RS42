function [error] = cal_annualerr(x, optV_tmp, model, method, spx_cell, vix_cell, upV, total_num)
% calculating the annual errors.
diff_list = cal_annualerr_lm(x, optV_tmp, model, method, spx_cell, vix_cell, upV, total_num);
error = sqrt(sum(diff_list.^2));
end

