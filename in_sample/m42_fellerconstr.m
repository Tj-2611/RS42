% This code estimates the H2F and 4/2SV model with the Feller condition as a constriant.
% Please first download "Differential Evolution" by Markus Buehren according to the url in "README.md"

clear;
currentDir = pwd;
parentDir = fileparts(currentDir);
addpath([parentDir, '/DE/'])
addpath([parentDir, '/source/'])
data_path = [parentDir, '/data/'];
setrandomseed(2024);
% model = '42withmu';
model = 'H2Fwithmu';
method = 'COS';
save_path = [parentDir, '/estimates/'];
if exist(save_path, 'dir') == 0
    mkdir(save_path)
end

parpool('local', 26);
retol = 1e-4;

year = 2016;
tic;
[optV_all, rmseV_all, optmdl_all, rmsemdl_all] = joca_1w1ysep_feller(year, data_path, save_path, model, method, retol);
disp(['Elapse time in year ', num2str(year), ' is ', num2str(toc), ' seconds'])

delete(gcp('nocreate'));

