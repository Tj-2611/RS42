% Out-of-sample test of the models.
% Please first run the corresponding in-sample codes to obtain the parameter estimates.

clear
currentDir = pwd;
parentDir = fileparts(currentDir);
addpath([parentDir, '/DE/'])
addpath([parentDir, '/source/'])
data_path = [parentDir, '/data/'];

% model = '42withmu';
% model = 'H2Fwithmu';
model = '42rswithmu';
method = 'COS';
save_path = 'workspace/';
if exist(save_path, 'dir') == 0
    mkdir(save_path)
end

tic;
parpool('local', 26);
hy = nan;
for year = 2014:2019
       [optV_all, rmseV_all, optmdl_all, rmsemdl_all] =os_joca_1w1ysep(year, hy, data_path, save_path, model, method);
end
 
year = 2020;
for hy = 1:2
       [optV_all, rmseV_all, optmdl_all, rmsemdl_all] = os_joca_1w1ysep(year, hy, data_path, save_path, model, method);
end
delete(gcp('nocreate'));
toc;
