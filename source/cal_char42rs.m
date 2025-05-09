function [phi] = cal_char42rs(params, tau, N, a, b, tgrid_num, expiries)
% calculating the characteristic functions under the 4/2RS model

k   = 0:N-1;
u   = k*pi/(b-a);

T = max(expiries);
V1 = params(1); 
H = params(4);
k1 = params(5);
theta1 = params(6);
sigma1 = params(7);
alpha = H + 0.5;  

dt = T/tgrid_num;
H = min(H, 0.4999);
[xn_ary, wn_ary] = markovianappr(H, T, dt, tgrid_num, 20); % Markovian approximation.

t_num = tgrid_num;
bb = -k1;

Tall_ary = 0:dt:(T+tau);
g0all_ary = V1 + k1*theta1/gamma(alpha+1)*Tall_ary.^alpha;

tau_num = length(Tall_ary(Tall_ary<tau));
E1all_ary = zeros(tau_num, 1);
for t_index = 1:tau_num
    s_tmp = Tall_ary(t_index);
    N_tmp = 1e3;
    n_ary = 0:N_tmp;
    ttmp = bb.^n_ary./gamma(alpha*(n_ary+1)+1).*s_tmp.^(alpha*(n_ary + 1));
    E1all_ary(t_index) = sum(ttmp);    
end

NN = length(xn_ary);
u_num = size(u, 2);
F_mt = zeros(t_num, u_num)*nan;
N_tmp = 1e3;
n_ary = 0:N_tmp;
initial_ary = bb.^n_ary.*tau.^(alpha*(n_ary+1)-1)./gamma(alpha*(n_ary+1)+1);
for u_index = 1:u_num
    phi_ary = zeros(t_num, NN);
    phi_ary(1, :) = 1i*u(u_index)*sum(initial_ary);
    for t_index = 2:t_num+1
        sum_tmp = wn_ary' * phi_ary(t_index-1, :).';
        F_mt(t_index-1, u_index) = -k1*sum_tmp + 1/2*sigma1^2*sum_tmp^2;
        phi_ary(t_index, :) = (phi_ary(t_index-1, :) + dt*F_mt(t_index-1))./(1+xn_ary'*dt);
    end
end


L2 = length(expiries);
phi = zeros(L2, u_num);
for i = 1:L2
    [~, ind_tmp] = min(abs(Tall_ary - expiries(i)));
    gT_ary = g0all_ary(ind_tmp:(ind_tmp+tau_num-1));
    integ1 = flip(gT_ary)*(E1all_ary*bb + 1)/tau*dt * 1i*u;
    gt_ary = g0all_ary(1:(ind_tmp-1));
    integ2 = flip(gt_ary)*F_mt(1:(ind_tmp-1), :)*dt;
    phi(i, :) = exp(integ1 + integ2);
end

end

