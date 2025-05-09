function [V_ary] = sim_rough(k, theta, sigma, V0, xn_ary, wn_ary, dt, t_num)
% simulate the rough component

n = length(xn_ary);
V_ary = zeros(t_num+1, 1);
V_ary(1) = V0;
Vi_ary = zeros(n, t_num+1);
noise = normrnd(0, 1, [t_num, 1]);
for t_index = 2:t_num+1
    Vplus = max(V_ary(t_index - 1), 0);
    Vi_ary(:, t_index) = exp(-xn_ary(:)*dt).*(Vi_ary(:, t_index-1) + ...
            k*(theta - Vplus)*dt + sigma*sqrt(Vplus)*noise(t_index-1)*sqrt(dt));
    V_ary(t_index) = wn_ary'*Vi_ary(:, t_index) + V0;
end

end