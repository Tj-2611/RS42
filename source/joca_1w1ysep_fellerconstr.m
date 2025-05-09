function [optV_all, rmseV_all, optmdl_all, rmsemdl_all] = joca_1w1ysep_fellerconstr(year, data_path, save_path, model, method, retol)
% parameter esitmations using two-step iterative algorithm 
% with the Feller condition as a constraint
filename = [num2str(year), '_', model, '_Fellerconstr.mat'];
savevars = {"optV_all", "rmseV_all", "optmdl_all", "rmsemdl_all", "year",  "model", "spxoption", "vixoption", "optV_tmp"};
opt = 'lm';

spxoption = readtable([data_path, 'surf_filt_liq.csv']);
spxoption = spxoption(spxoption.year == year, :);
vixoption = readtable([data_path, 'vix_option_future.csv']); 
vixoption = vixoption(vixoption.year == year, :);

spxoption.log_mon = log(spxoption.K./(spxoption.S.*exp((spxoption.interest_rate - spxoption.dividend_rate).*spxoption.t)));
spxoption = spxoption(spxoption.log_mon > 0, :);
spxoption = spxoption(spxoption.log_mon < log(1.25), :);
spxoption = spxoption(spxoption.t >= 7/365, :);
spxoption = spxoption(spxoption.t <= 1, :);
spxoption = spxoption(spxoption.volume >= 5, :);
vixoption = vixoption(vixoption.log_mon > 0, :);
vixoption = vixoption(vixoption.log_mon <  log(1.25), :);
vixoption = vixoption(vixoption.t >= 7/365, :);
vixoption = vixoption(vixoption.t <= 1, :);
vixoption = vixoption(vixoption.volume >= 5, :);

spxoption = spxoption(spxoption.best_bid>3/8, :);
vixoption = vixoption(vixoption.best_bid>3/8, :);

spx_variables_cell = {'date', 't', 'K', 'mid_price', 'best_bid', 'best_offer', 'impl_volatility', 'interest_rate', 'dividend_rate', 'S', 'vega'};
vix_variables_cell = {'date', 't', 'K', 'mid_price', 'best_bid', 'best_offer', 'impl_volatility', 'interest_rate', 'S', 'vega'};
spxoption = spxoption(:, spx_variables_cell);
vixoption = vixoption(:, vix_variables_cell);


date_list = unique([spxoption.date; vixoption.date]);
t_num = size(date_list, 1);
M = size(spxoption, 1);
MM = size(vixoption, 1);
spx_cell = cell(t_num, 1); 
vix_cell = cell(t_num, 1);
for i = 1:t_num
    spx_cell{i} = spxoption(spxoption.date == date_list(i), :);
    vix_cell{i} = vixoption(vixoption.date == date_list(i), :);
end

[paramdl0, paramdl_bd, options, optimInfo, paramDefCell, ...
    objFctSettings, objFctParams, DEParams, emailParams] = initial_opt(model, opt);
[V0,V_bd, options_V, optimInfo_V, paramDefCell_V, ...
    objFctSettings_V, objFctParams_V, DEParams_V, emailParams_V] = initial_V(model, opt);

V_num = size(V0, 1);
para_num = size(paramdl0, 1);
step_num = 1e2;
optV_all = zeros(V_num, t_num, step_num);
rmseV_all = zeros(t_num, step_num);
optmdl_all = zeros(para_num, step_num);
rmsemdl_all = zeros(step_num, 1);

optV_tmp = zeros(V_num, t_num)*nan;
for i = 1:V_num
    optV_tmp(i, :) = V0(i);
end
paramdl_tmp = paramdl0;


tol_tmp = 1;
iter = 1;
% 
while tol_tmp > retol
    tic;
    rmseV_tmp = zeros(t_num, 1)*nan;
   if iter ~=1
   for t_index = 1:t_num
      V_tmp = optV_tmp(:, t_index);
      spx_tmp = spx_cell{t_index};
      vix_tmp = vix_cell{t_index};
      if iter == 2
          paramDefCell_V = {'', V_bd, zeros(length(V0), 1), V_tmp};  
          errV_de = @(x)(cal_volsurfsep([x; paramdl_tmp], spx_tmp, vix_tmp, model, method, 'ip'));
          [bestmem, ~, ~, ~, ~] = differentialevolution(...
                    DEParams_V, paramDefCell_V, errV_de, objFctSettings_V, objFctParams_V, emailParams_V, optimInfo_V); 
          V_tmp = bestmem;
      end
      errV = @(x)(cal_volsurfsep([x; paramdl_tmp], spx_tmp, vix_tmp, model, method, 'lm'));
      [V_tmp,fval] = lsqnonlin(errV,V_tmp,V_bd(:, 1),V_bd(:, 2),options);
      rmse = sqrt(fval/(size(spx_tmp, 1)+size(vix_tmp, 1)));
      optV_tmp(:, t_index) = V_tmp;
      rmseV_tmp(t_index) = rmse;
      fprintf('<<OptV: The optimized RMSE is %f\n', rmse)
   end
    optV_all(:, :, iter) = optV_tmp;
    rmseV_all(:, iter) = rmseV_tmp;
   end
   
    errmdl = @(x)(cal_annualerr_lm(x, optV_tmp, model, method, spx_cell, vix_cell, M+MM));
    if strcmp(model, '42withmu')
       [paramdl_tmp, fval] = lsqnonlin(errmdl,paramdl_tmp,paramdl_bd(:, 1),paramdl_bd(:, 2),[],[],[],[],Feller_constraint42,options);
    elseif  strcmp(model, 'H2Fwithmu')
       [paramdl_tmp, fval] = lsqnonlin(errmdl,paramdl_tmp,paramdl_bd(:, 1),paramdl_bd(:, 2),[],[],[],[],Feller_constraintH2F,options);  
    else
       [paramdl_tmp, fval] = lsqnonlin(errmdl,paramdl_tmp,paramdl_bd(:, 1),paramdl_bd(:, 2),options);          
    end    
    rmse1 = sqrt(fval);
    optmdl_all(:, iter) = paramdl_tmp;
    rmsemdl_all(iter, 1) = rmse1;
    
    save([save_path, '/', filename], savevars{:});
    if iter > 1
        tol_tmp = abs(rmsemdl_all(iter)/rmsemdl_all(iter-1) - 1);
        fprintf('<<year: %d, iteration: %d/%d, optimized RMSE is %.5f, relative tolerance is %.5f\n', year, iter, step_num, rmse1, tol_tmp)  
    else
        fprintf('<<year: %d, iteration: %d/%d, optimized RMSE is %.5f\n', year, iter, step_num, rmse1)
    end
    iter = iter + 1;
    if iter > 100
        break;
    end
end
save([save_path, '/', filename], savevars{:});
end


function [c, ceq] = Feller_constraint42(x)
k1 = x(2);
theta1 = x(3);
sigma1 = x(4);
c = sigma1^2 - 2*k1*theta1;
ceq = [];
end



function [c, ceq] = Feller_constraintH2F(x)
k1 = x(2);
theta1 = x(3);
sigma1 = x(4);

k2 = x(6);
theta2 = x(7);
sigma2 = x(8);

c1 = sigma1^2 - 2*k1*theta1;
c2 = sigma2^2 - 2*k2*theta2;
c = [c1; c2];
ceq = [];
end

