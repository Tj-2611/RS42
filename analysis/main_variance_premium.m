% Calculating variance_risk_premium and visualisation.
clear;
currentDir = pwd;
parentDir = fileparts(currentDir);
addpath([parentDir, '/source/'])
filepath = [parentDir, '/estimates/'];
data_path = [parentDir, '/data/'];
year_list = [2013:2019, 2020:0.5:2020.5];
analysis_model = '42rswithmu';
rsV = [];
for year = year_list
    if year >= 2020
        yr = floor(year);
        hy = ((year-yr)/0.5 + 1);
    else
        yr = year;
        hy = [];
    end
    [model, filename] = find_model(yr, hy, analysis_model); 
    load([filepath, filename], 'optV_all', 'optmdl_all', 'spxoption', 'vixoption', 'rmseV_all');
    iter = sum(any(optmdl_all~=0, 1));
    optV = optV_all(:, :, iter);
    paramdl = optmdl_all(:, iter);

    date_list = unique([spxoption.date; vixoption.date]);
    rsV_tb = table();
    rsV_tb.date = date_list;
    rsV_tb.rmse = rmseV_all(:, iter);
    rsV_tb.V1 = optV(1, :)';
    rsV_tb.V2 = optV(2, :)';
    mu = paramdl(1);
    rsV_tb.V = mu*rsV_tb.V1 + (1-mu)*rsV_tb.V2;
    rsV = [rsV; rsV_tb];
%    
end
rv_tb = readtable([data_path, 'RV.csv']);
rv_tb = table2timetable(rv_tb,'RowTimes','Date');
realVar = [];
vix_tb =  readtable([data_path, 'cboe_vix.csv']);
vix_tb = table2timetable(vix_tb,'RowTimes','Date');
Vix = [];
window = 22;

for i = 1:length(rsV.date)
    dt = rsV.date(i);
    tmp = rv_tb(rv_tb.Date <= dt, :);
    tmp = tmp(strcmp(tmp.Type, 'QMLE-Trade'), 'Volatility');
    tmp = tmp(end-window+1: end, :);
    realVar_tmp = sum(tmp.Volatility.^2)/22;
    ttmp = vix_tb(vix_tb.Date <= dt, :);
    Vix = [Vix; ttmp.vix(end)];
    realVar = [realVar; realVar_tmp];
end

realVix_tb = table();
realVix_tb.realVar = realVar;
realVix_tb.vix = Vix;
realVix_tb.date = rsV.date;
realVix_tb = table2timetable(realVix_tb,'RowTimes','date');

y = realVix_tb.realVar(2:end);
x1 = realVix_tb.realVar(1:end-1);
x2 = (realVix_tb.vix(1:end-1)/100).^2;

X = [x1, x2];
mdl = fitlm(X, y);

expVar = predict(mdl,X);
realVix_tb.expVar = [nan; expVar];

du = 5e-2;
rsV = [];
for year = year_list
    if year >= 2020
        yr = floor(year);
        hy = ((year-yr)/0.5 + 1);
    else
        yr = year;
        hy = [];
    end
    [model, filename] = find_model(yr, hy, analysis_model); 
    load([filepath, filename], 'optV_all', 'optmdl_all', 'spxoption', 'vixoption', 'rmseV_all');
    iter = sum(any(optmdl_all~=0, 1));
    optV = optV_all(:, :, iter);
    paramdl = optmdl_all(:, iter);
    date_list = unique([spxoption.date; vixoption.date]);
    EQV_ary = zeros(length(date_list),1);
    for  ii = 1:length(date_list)
         EQV_tmp = expQv(model, [optV(:, ii); paramdl], du);
         EQV_ary(ii) = EQV_tmp;
    end
    rsV_tb = table();
    rsV_tb.date = date_list;
    rsV_tb.rmse = rmseV_all(:, iter);
    rsV_tb.V1 = optV(1, :)';
    rsV_tb.V2 = optV(2, :)';
    rsV_tb.EQV = EQV_ary;

    mu = paramdl(1);
    rsV_tb.V = mu*rsV_tb.V1 + (1-mu)*rsV_tb.V2;
    rsV = [rsV; rsV_tb];
end


nonaffV = [];
analysis_model = '42withmu';
for year = year_list
    if year >= 2020
        yr = floor(year);
        hy = ((year-yr)/0.5 + 1);
    else
        yr = year;
        hy = [];
    end
    [model, filename] = find_model(yr, hy, analysis_model); 
    load([filepath, filename], 'optV_all', 'optmdl_all', 'spxoption', 'vixoption', 'rmseV_all');
    iter = sum(any(optmdl_all~=0, 1));
    optV = optV_all(:, :, iter);
    paramdl = optmdl_all(:, iter);
    date_list = unique([spxoption.date; vixoption.date]);
    EQV_ary = zeros(length(date_list),1);
    for  ii = 1:length(date_list)
         EQV_tmp = expQv(model, [optV(:, ii); paramdl], du);
         EQV_ary(ii) = EQV_tmp;
    end
    rsV_tb = table();
    rsV_tb.date = date_list;
    rsV_tb.rmse = rmseV_all(:, iter);
    rsV_tb.V1 = optV(1, :)';
    rsV_tb.V2 = optV(2, :)';
    rsV_tb.EQV = EQV_ary;
    mu = paramdl(1);
    rsV_tb.V = mu*rsV_tb.V1 + (1-mu)*rsV_tb.V2;
    nonaffV = [nonaffV; rsV_tb];  
end


affV = [];
analysis_model = 'H2Fwithmu';
for year = year_list
    if year >= 2020
        yr = floor(year);
        hy = ((year-yr)/0.5 + 1);
    else
        yr = year;
        hy = [];
    end
    [model, filename] = find_model(yr, hy, analysis_model); 
    load([filepath, 'H2Fsep/', filename], 'optV_all', 'optmdl_all', 'spxoption', 'vixoption', 'rmseV_all');
    iter = sum(any(optmdl_all~=0, 1));
    optV = optV_all(:, :, iter);
    paramdl = optmdl_all(:, iter);
    date_list = unique([spxoption.date; vixoption.date]);
    EQV_ary = zeros(length(date_list),1);
    for  ii = 1:length(date_list)
         EQV_tmp = expQv(model, [optV(:, ii); paramdl], du);
         EQV_ary(ii) = EQV_tmp;
    end
    rsV_tb = table();
    rsV_tb.date = date_list;
    rsV_tb.rmse = rmseV_all(:, iter);
    rsV_tb.V1 = optV(1, :)';
    rsV_tb.V2 = optV(2, :)';
    rsV_tb.EQV = EQV_ary;
    mu = paramdl(1);
    rsV_tb.V = mu*rsV_tb.V1 + (1-mu)*rsV_tb.V2;
    affV = [affV; rsV_tb];
%    
end

rsV.VP = realVix_tb.expVar - rsV.EQV;
nonaffV.VP = realVix_tb.expVar- nonaffV.EQV;
affV.VP = realVix_tb.expVar - affV.EQV;


%%
mfV = table();
mfV.date = realVix_tb.date;
mfV.EQV = (realVix_tb.vix/100).^2;


%%
figure;
plot(rsV.date(2:end), realVix_tb.expVar(2:end) - rsV.EQV(2:end))
hold on;
plot(nonaffV.date(2:end), realVix_tb.expVar(2:end) - nonaffV.EQV(2:end), 'LineStyle','-.')
hold on;
plot(affV.date(2:end), realVix_tb.expVar(2:end) - affV.EQV(2:end), 'LineStyle','--')
hold on;
plot(mfV.date(2:end), realVix_tb.expVar(2:end) - mfV.EQV(2:end), 'LineStyle',':', color = [0.5, 0.5, 0.5])
legend({"4/2RS", "4/2SV", "H2F", "Benchmark"}, 'FontSize',13)
xlabel('Time', 'FontSize',13)
ylabel('Variance risk premium', 'FontSize',13)
datetick('x', 'yyyy');


mean(realVix_tb.expVar(2:end) - rsV.EQV(2:end))
% -0.0121
mean(realVix_tb.expVar(2:end) - nonaffV.EQV(2:end))
% -0.0111
mean(realVix_tb.expVar(2:end) - affV.EQV(2:end))
% -0.0261
mean(realVix_tb.expVar(2:end) - mfV.EQV(2:end))
% -0.0182

