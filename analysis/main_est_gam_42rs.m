% Estimating EVs using the method by Aıt-Sahalia and Kimmel (2007)
% 4/2RS model
clear;
addpath('ak07method');
rng('default');
save_path = 'workspace\';
if exist(save_path, 'dir') == 0
    mkdir(save_path)
end


mu_list = linspace(0, 0.95, 10);
mu_num = length(mu_list);

paramdl_tmp = [0.2; 0.2; 1; 0.5; 1;-0.7;
                         50; 1; 50; -0.7];       
params = [0.1; 0.1; paramdl_tmp];                                       
n = 1.5e3;
tic;
del = 1/252;
factor = 100;
nsimul = factor*n;
delsimul = del/factor;
burnin = 5e2;

V1 = params(1);
V2 = params(2);
mu = params(3);
H = params(4);
k1 = params(5);
theta1 = params(6);
sigma1 = params(7);
k2 = params(9);
theta2 = params(10);
sigma2 = params(11);
T = del*n;
[xn_ary, wn_ary] = markovianappr(H, T, delsimul, nsimul, 15);

parpool('local', 26);
rep_num = 1e3;
rng('default') 
astart=0.1; bstart=0.1;  cstart=50; dstart = 1.5;
est_mt = zeros(4, rep_num, mu_num);
se_mt = zeros(4, rep_num, mu_num);
flag_mt = zeros(rep_num, mu_num);
for ii = 1:mu_num
    mu_tmp = mu_list(ii);
    param0 = [astart,bstart,cstart,dstart]';
    parfor rep = 1:rep_num
        flag = 0;
        [V1_ary] = sim_rough(k1, theta1, sigma1, V1, xn_ary, wn_ary, delsimul, nsimul);
        [V2_ary] = sim_32(k2, theta2, sigma2, V2, delsimul, nsimul);
        V_ary = mu_tmp*V1_ary + (1-mu_tmp)*V2_ary;
        V_ary = V_ary((burnin*factor + 1):end);
        x = zeros(n-burnin,1);
        for i = 1:n-burnin
            x(i) = V_ary(1+(i-1)*factor);
        end
        output = mymle(@ModelU3, x, del, param0);

        est_mt(:, rep, ii) = output.param;
        se_mt(:, rep, ii) = output.se;
        flag_mt(rep, ii) = output.exitflag;
    end
    fprintf('<< step %d/%d \n', ii, mu_num);
end
delete(gcp('nocreate'));

save('workspace/varymu_1e3_42rs.mat')
%%
% figure;
% est_mt_tmp = est_mt;
% gam_mt = reshape(est_mt_tmp(4, :, :), rep_num, mu_num);
% bp = boxplot(gam_mt, 'Symbol', '');
% xticklabels(round(mu_list, 2))    
% xlabel('Proportion of the component V_1', 'FontSize', 13)
% ylabel('Elasticity of variance',  'FontSize', 13)

