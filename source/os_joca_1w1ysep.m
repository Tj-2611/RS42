function [optV_all, rmseV_all, optmdl_all, rmsemdl_all] = os_joca_1w1ysep(year, hy, data_path, save_path, analysis_model, method)
% out-of-sample test by calibrating spot variances with fixed structural parameters 
currentDir = pwd;
parentDir = fileparts(currentDir);
parafilepath = [parentDir, '/estimates/'];
if year == 2021
    [model, parafilename] = find_model(year-1, 2, analysis_model);
elseif year == 2020
    if hy == 2
        [model, parafilename] = find_model(year, 1, analysis_model);
    else
        [model, parafilename] = find_model(year-1, 1, analysis_model);
    end
else
    [model, parafilename] = find_model(year-1, hy, analysis_model);
end
load([parafilepath, parafilename],'optmdl_all');
iter = sum(any(optmdl_all~=0, 1));
paramdl = optmdl_all(:, iter);            
            
opt = 'lm';
spxoption = readtable([data_path, 'surf_filt_liq.csv']);
spxoption = spxoption(spxoption.year == year, :);
vixoption = readtable([data_path, 'vix_option_future.csv']); 
vixoption = vixoption(vixoption.year == year, :);

if year == 2021 || year == 2020
    if hy == 1
        spxoption.month = month(spxoption.date);
        spxoption = spxoption(spxoption.month <= 6, :); 
        vixoption.month = month(vixoption.date);
        vixoption = vixoption(vixoption.month <= 6, :);
    else
        spxoption.month = month(spxoption.date);
        spxoption = spxoption(spxoption.month > 6, :); 
        vixoption.month = month(vixoption.date);
        vixoption = vixoption(vixoption.month > 6, :);
    end
else
    spxoption = spxoption(spxoption.year == year, :);
    vixoption = vixoption(vixoption.year == year, :);
end

    
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

options_tmp1 = optimoptions('fmincon','Algorithm','interior-point',...
'MaxFunctionEvaluations',10^4,'MaxIterations',30,'FunctionTolerance',10^(-12),...
'StepTolerance',10^(-4),'FiniteDifferenceStepSize',10^(-4),'Display','iter-detailed');
options_tmp2 = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt',...
'MaxFunctionEvaluations',10^4,'MaxIterations',30,'FunctionTolerance',10^(-12),...
'StepTolerance',10^(-4),'FiniteDifferenceStepSize',10^(-4),'Display','iter-detailed');

[V0, V_bd, ~, ~, ~, ~, ~, ~, ~] = initial_V(model, opt);

V_num = size(V0, 1);
step_num = 1;
para_num = size(paramdl, 1);
optmdl_all = zeros(para_num, step_num);


optV_tmp = zeros(V_num, t_num)*nan;
for i = 1:V_num
    optV_tmp(i, :) = V0(i);
end

paramdl_tmp = paramdl;

%%
rmseV_tmp = zeros(t_num, 1)*nan;
squareerro_list = zeros(t_num, 1);
parfor t_index = 1:t_num
  V_tmp = optV_tmp(:, t_index);
  spx_tmp = spx_cell{t_index};
  vix_tmp = vix_cell{t_index};
  errV1 = @(x)(cal_volsurfsep([x; paramdl_tmp], spx_tmp, vix_tmp, model, method, 'ip'));
  [V_tmp1,fval1] = fmincon(errV1,[V0],[],[],[],[],[V_bd(:, 1)],[V_bd(:, 2)],[], options_tmp1);
  errV2 = @(x)(cal_volsurfsep([x; paramdl_tmp], spx_tmp, vix_tmp, model, method, 'lm'));
  [V_tmp2,fval2] = lsqnonlin(errV2,[V0],[V_bd(:, 1)],[V_bd(:, 2)],options_tmp2);  
  if fval1 >= fval2 
      fval = fval2;
      V_tmp = V_tmp2;
  else
      fval = fval1;
      V_tmp = V_tmp1;     
  end
  rmse = sqrt(fval/(size(spx_tmp, 1)+size(vix_tmp, 1)));
  squareerro_list(t_index) = fval
  optV_tmp(:, t_index) = V_tmp;
  rmseV_tmp(t_index) = rmse;
  fprintf('<<OptV: The optimized RMSE is %f\n', rmse)
end
optV_all = optV_tmp;
rmseV_all = rmseV_tmp;
rmsemdl_all = sqrt(sum(squareerro_list)/(M+MM));

if year < 2020
    filename = ['os', num2str(year), '_', analysis_model, '.mat'];
elseif  year >= 2020
    filename = ['os', num2str(year), '_half', num2str(hy), '_', analysis_model, '.mat'];
end

savevars = {"optV_all", "rmseV_all", "paramdl", "rmsemdl_all",...
            "year",  "model", "spxoption", "vixoption"};

fprintf('<<year: %d, hy: %d, optimized RMSE is %.5f\n', year, hy, rmsemdl_all)
save([save_path, '/', filename], savevars{:});
end


