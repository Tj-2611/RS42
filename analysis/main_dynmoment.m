% Calculating risk-neutrual moments of log-return
clear;

currentDir = pwd;
parentDir = fileparts(currentDir);
addpath([parentDir, '/source/'])
filepath = [parentDir, '/estimates/'];
data_path = [parentDir, '/data/'];
spxoption = readtable([data_path, 'surf_filt_liq.csv']);
vixoption = readtable([data_path, 'vix_option_future.csv']); 
spxoption.log_mon = log(spxoption.K./(spxoption.S.*exp((spxoption.interest_rate - spxoption.dividend_rate).*spxoption.t)));
spxoption = spxoption(spxoption.log_mon > 0, :);
spxoption = spxoption(spxoption.log_mon < log(1.25), :);
spxoption = spxoption(spxoption.t >= 7/365, :);
spxoption = spxoption(spxoption.t <= 1, :);
spxoption = spxoption(spxoption.volume >= 5, :);
vixoption = vixoption(vixoption.log_mon > 0, :);
vixoption = vixoption(vixoption.log_mon < log(1.25), :);
vixoption = vixoption(vixoption.t >= 7/365, :);
vixoption = vixoption(vixoption.t <= 1, :);
vixoption = vixoption(vixoption.volume >= 5, :);

spxoption = spxoption(spxoption.best_bid>3/8, :);
vixoption = vixoption(vixoption.best_bid>3/8, :);

year_list = [2013:2019, 2020:0.5:2020.5];

analysis_model = '42rswithmu';
T = 1/12;
du = 1e-2;
ex_list = [];
var_list = [];
skew_list = [];
kurt_list = [];
date_list = [];
all_listtmp = [];
for year = year_list
    if year >= 2020
        yr = floor(year);
        hy = ((year-yr)/0.5 + 1);
    else
        yr = year;
        hy = [];
    end
    [model, filename] = find_model(yr, hy, analysis_model); 
    load([filepath, filename], 'optV_all', 'optmdl_all', 'spxoption', 'vixoption');
    
    iter = sum(any(optmdl_all~=0, 1));
    optV = optV_all(:, :, iter);
    paramdl = optmdl_all(:, iter);
    dt = unique([spxoption.date; vixoption.date]);
    date_list = [date_list; dt];
    for j = 1:size(optV, 2)
        params = [optV(:, j); paramdl];
        dt_tmp = dt(j);
        IRtb = spxoption(spxoption.date==dt_tmp, 'interest_rate');
        IR = IRtb.interest_rate(1);
        [Ex, Var, Skew, Kurt] = cal_4moment(params, model, du, IR, T);
        all_listtmp = [all_listtmp; [Ex, Var, Skew, Kurt]];
        ex_list = [ex_list; Ex];
        var_list = [var_list; Var];
        skew_list = [skew_list; Skew];
        kurt_list = [kurt_list; Kurt];
    end
end

mors_tb = table();
mors_tb.date = date_list;
mors_tb.ex = ex_list;
mors_tb.var = var_list;
mors_tb.skew = skew_list;
mors_tb.kurt = kurt_list;
mors_tb = table2timetable(mors_tb, "RowTimes","date");


analysis_model = '42withmu';
T = 1/12;
du = 1e-2;
ex_list = [];
var_list = [];
skew_list = [];
kurt_list = [];
date_list = [];
all_listtmp = [];
for year = year_list
    if year >= 2020
        yr = floor(year);
        hy = ((year-yr)/0.5 + 1);
    else
        yr = year;
        hy = [];
    end
    [model, filename] = find_model(yr, hy, analysis_model); 
    load([filepath, filename], 'optV_all', 'optmdl_all', 'spxoption', 'vixoption');
    iter = sum(any(optmdl_all~=0, 1));
    optV = optV_all(:, :, iter);
    paramdl = optmdl_all(:, iter);
    dt = unique([spxoption.date; vixoption.date]);
    date_list = [date_list; dt];
    for j = 1:size(optV, 2)
        params = [optV(:, j); paramdl];
        dt_tmp = dt(j);
        IRtb = spxoption(spxoption.date==dt_tmp, 'interest_rate');
        IR = IRtb.interest_rate(1);
        [Ex, Var, Skew, Kurt] = cal_4moment(params, model, du, IR, T);
        all_listtmp = [all_listtmp; [Ex, Var, Skew, Kurt]];
        ex_list = [ex_list; Ex];
        var_list = [var_list; Var];
        skew_list = [skew_list; Skew];
        kurt_list = [kurt_list; Kurt];
    end
end

mo42_tb = table();
mo42_tb.date = date_list;
mo42_tb.ex = ex_list;
mo42_tb.var = var_list;
mo42_tb.skew = skew_list;
mo42_tb.kurt = kurt_list;
mo42_tb = table2timetable(mo42_tb, "RowTimes","date");



analysis_model = 'H2Fwithmu';
T = 1/12;
du = 1e-2;
ex_list = [];
var_list = [];
skew_list = [];
kurt_list = [];
date_list = [];
all_listtmp = [];
for year = year_list
    if year >= 2020
        yr = floor(year);
        hy = ((year-yr)/0.5 + 1);
    else
        yr = year;
        hy = [];
    end
    [model, filename] = find_model(yr, hy, analysis_model); 
    load([filepath, 'H2Fsep\', filename], 'optV_all', 'optmdl_all', 'spxoption', 'vixoption');
    iter = sum(any(optmdl_all~=0, 1));
    optV = optV_all(:, :, iter);
    paramdl = optmdl_all(:, iter);
    dt = unique([spxoption.date; vixoption.date]);
    date_list = [date_list; dt];
    for j = 1:size(optV, 2)
        params = [optV(:, j); paramdl];
        dt_tmp = dt(j);
        IRtb = spxoption(spxoption.date==dt_tmp, 'interest_rate');
        IR = IRtb.interest_rate(1);
        [Ex, Var, Skew, Kurt] = cal_4moment(params, model, du, IR, T);
        all_listtmp = [all_listtmp; [Ex, Var, Skew, Kurt]];
        ex_list = [ex_list; Ex];
        var_list = [var_list; Var];
        skew_list = [skew_list; Skew];
        kurt_list = [kurt_list; Kurt];
    end
end

moH2F_tb = table();
moH2F_tb.date = date_list;
moH2F_tb.ex = ex_list;
moH2F_tb.var = var_list;
moH2F_tb.skew = skew_list;
moH2F_tb.kurt = kurt_list;
moH2F_tb = table2timetable(moH2F_tb, "RowTimes","date");

datetick('x', 'yyyy-MM-dd');
figure;
plot(mors_tb.date, mors_tb.ex, '-')
hold on;
plot(mo42_tb.date, mo42_tb.ex, '-.')
hold on;
plot(moH2F_tb.date, moH2F_tb.ex, '--')
legend({"4/2RS", "4/2SV", "H2F"},  'FontSize', 13)
xlabel('Time', 'FontSize', 13)
ylabel('Mean', 'FontSize', 13)
datetick('x', 'yyyy');

figure;
plot(mors_tb.date, mors_tb.var, '-')
hold on;
plot(mo42_tb.date, mo42_tb.var, '-.')
hold on;
plot(moH2F_tb.date, moH2F_tb.var, '--')
hold on;
legend({"4/2RS", "4/2SV", "H2F"},  'FontSize', 13)
xlabel('Time',  'FontSize', 13)
ylabel('Variance',  'FontSize', 13)
datetick('x', 'yyyy');

% 
figure;
plot(mors_tb.date, mors_tb.skew, '-')
hold on;
plot(mo42_tb.date, mo42_tb.skew, '-.')
hold on;
plot(moH2F_tb.date, moH2F_tb.skew, '--')
legend({"4/2RS", "4/2SV", "H2F"},  'FontSize', 13)
xlabel('Time',  'FontSize', 13)
ylabel('Skewness',  'FontSize', 13)
datetick('x', 'yyyy');

figure;
plot(mors_tb.date, mors_tb.kurt, '-')
hold on;
plot(mo42_tb.date, mo42_tb.kurt, '-.')
hold on;
plot(moH2F_tb.date, moH2F_tb.kurt, '--')
legend({"4/2RS", "4/2SV", "H2F"},  'FontSize', 13)
xlabel('Time',  'FontSize', 13)
ylabel('Excess kurtosis',  'FontSize', 13)
datetick('x', 'yyyy');



figure;
subplot(2,2,1)
plot(mo42_tb.date, mo42_tb.ex)
hold on;
plot(mors_tb.date, mors_tb.ex)
hold on;
plot(moH2F_tb.date, moH2F_tb.ex)
subplot(2,2,2)
plot(mo42_tb.date, mo42_tb.var)
hold on;
plot(mors_tb.date, mors_tb.var)
hold on;
plot(moH2F_tb.date, moH2F_tb.var)
subplot(2,2,3)
plot(mo42_tb.date, mo42_tb.skew)
hold on;
plot(mors_tb.date, mors_tb.skew)
hold on;
plot(moH2F_tb.date, moH2F_tb.skew)
subplot(2,2,4)
plot(mo42_tb.date, mo42_tb.kurt)
hold on;
plot(mors_tb.date, mors_tb.kurt)
hold on;
plot(moH2F_tb.date, moH2F_tb.kurt)
legend({"4/2SV", "4/2RS", "H2F"})




