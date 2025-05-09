% This code estimates the H2F model based on SPX options only.
% Please download "Differential Evolution" by Markus Buehren according to the url in "README.md"

clear
currentDir = pwd;
parentDir = fileparts(currentDir);
addpath([parentDir, '/DE/']);
addpath([parentDir, '/source/']);
data_path = [parentDir, '/data/'];
setrandomseed(2024);

model = 'H2Fwithmu';
method = 'COS';
save_path = [parentDir, '/estimates/H2FwithmuSPXonly'];
if exist(save_path, 'dir') == 0
    mkdir(save_path)
end

parpool('local', 12);
retol = 1e-4;
for year = 2013:2019
	tic; 
    [optV_all, rmseV_all, optmdl_all, rmsemdl_all] = joca_1w1ysepSPXonly(year, data_path, save_path, model, method, retol);
	disp(['Elapse time in year ', num2str(year), 'num2str(hy)', ' is ', num2str(toc), ' seconds'])
end

year = 2020;
for hy = 1:2
% hy represents the first/second half of the year
    tic; 
    [optV_all, rmseV_all, optmdl_all, rmsemdl_all] = joca_1w1ysephySPXonly(year, hy, data_path, save_path, model, method, retol);
    disp(['Elapse time in year ', num2str(year), 'num2str(hy)', ' is ', num2str(toc), ' seconds'])
end

delete(gcp('nocreate'));
