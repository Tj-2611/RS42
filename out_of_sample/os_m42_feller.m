% Out-of-sample test of the H2F and 4/2SV model with the Feller constraint.
% Please first run the corresponding in-sample codes to obtain the parameter estimates.
clear
currentDir = pwd;
parentDir = fileparts(currentDir);
addpath([parentDir, '/DE/'])
addpath([parentDir, '/source/'])
data_path = [parentDir, '/data/'];

% model = 'H2Fwithmu';
model = '42withmu';
method = 'COS';
save_path = 'workspace/';
if exist(save_path, 'dir') == 0
    mkdir(save_path)
end

tic;
parpool('local', 26);

hy = nan;
for year = 2017
       [optV_all, rmseV_all, optmdl_all, rmsemdl_all] =os_joca_1w1ysep_fellerconstr(year, hy, data_path, save_path, model, method);
end

delete(gcp('nocreate'));
toc;
