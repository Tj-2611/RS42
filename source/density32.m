function [f2_ary] = density32(y_ary, T, k, sigma, theta, V0)
% transition density function of 3/2 component
delta = 4*(k+sigma^2)/sigma^2;

cT = sigma^2*(exp(k*theta*T)-1)/(4*k*theta);
alpha = 1/(V0*cT);

tily = exp(k*theta*T)./(y_ary*cT);
f2_ary = 1./y_ary.^2*exp(k*theta*T)/cT.*ncx2pdf(tily, delta, alpha);
f2_ary(1) = 0;
end

