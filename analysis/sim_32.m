function [V_ary] = sim_32(k, theta, sigma, V0, dt, t_num)
% simulate the 3/2 component
V_ary = zeros(t_num+1, 1);
V_ary(1) = V0;

delta = 4*(k+sigma^2)/sigma^2;
cdeltat = sigma^2*(exp(k*theta*dt)-1)/(4*k*theta);


for i = 2:(t_num+1)
    alphat = 1/(V_ary(i-1)*cdeltat);
    tmp = ncx2rnd(delta,alphat,1);
    V_ary(i) = exp(k*theta*dt)/cdeltat/tmp;
end


end


% backup
% V_ary = zeros(t_num+1, 1);
% V_ary(1) = V0;
% delta = 4*(k+sigma^2)/sigma^2;
% for i = 2:(t_num+1)
%     t = dt*(i-1);
%     ct = sigma^2*(exp(k*theta*t)-1)/(4*k*theta);
%     alphat = 1/(V0*ct);
%     tmp = ncx2rnd(delta,alphat,1);
%     V_ary(i) = 1./(tmp*ct/exp(k*theta*t));
% end


% y_num = 1e3;
% h0 = [0; 1]; 
% y_ary = linspace(1e-4, 0.5, y_num);
% [~,h_ary] = ode45(@(y,h) fun_h(y, h, k, sigma), y_ary*(exp(k*sigma*tau)-1)/(k*theta), h0);
% [f2_ary] = density32(y_ary, T, k, sigma, theta, V0);
% h_ary = h_ary(:, 1);