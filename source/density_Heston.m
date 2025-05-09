function [f2_ary] = density_Heston(x_ary, T, k, sigma, theta, V0)
% transition density of the Heston component
delta = 4*(k*theta)/sigma^2;

cT = sigma^2*(exp(k*T)-1)/(4*k);
alpha = V0/cT;

tilx = x_ary* exp(k*T)./cT;
f2_ary = exp(k*T)./cT.*ncx2pdf(tilx, delta, alpha);
f2_ary(1) = 0;
f2_ary = f2_ary';
end

