function [V_ary] = sim_12(k, theta, sigma, V0, dt, t_num)
% simulate the 1/2 component
V_ary = zeros(t_num+1, 1);
V_ary(1) = V0;

delta = 4*k*theta/sigma^2;
cdeltat = sigma^2*(exp(k*dt)-1)/(4*k);


for i = 2:(t_num+1)
    alphat = V_ary(i-1)/cdeltat;
    tmp = ncx2rnd(delta,alphat,1);
    V_ary(i) = tmp *cdeltat/exp(k*dt);
end


end

