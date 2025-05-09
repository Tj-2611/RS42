function [x_ary, w_ary] = markovianappr(H, T, dt, grid_num, NN)
% geometric Gaussian approximation method for approximating
% the fractional kernal
t_ary = linspace(dt, T, grid_num);
[x_mt, w_mt, ~] = cal_KN(NN, T, H, t_ary);
error_c = @(c) scale_K(c, NN, T, H, t_ary); 
[res, ~] = fminunc(error_c, 1, optimoptions(@fminunc,'Display','off'));
w_mt = w_mt*res;
x_mt = x_mt';
x_ary = x_mt(:);
w_mt = w_mt';
w_ary = w_mt(:);
end

function [L1norm] = scale_K(c, N, T, H, t_ary)

[~, ~, KNt] = cal_KN(N, T, H, t_ary);
K2_ary = t_ary.^(H-0.5)/gamma(H+0.5);
L1norm = mean(abs(K2_ary'-KNt*c))*(t_ary(end) - t_ary(1));
end


function [x_mt, w_mt, KNt] = cal_KN(N, T, H, t_ary)
aalpha = log(3+2*sqrt(2));
bbeta = 1;
aa = (10*sqrt(2)-14)/exp(1)*sqrt((H+0.5)*N)/T;
bb = (10*sqrt(2)-14)/exp(1)/T;

mm = round(bbeta*sqrt((0.5 +H)*N));
nn = round(1/bbeta*sqrt(N/(0.5 + H)));

xi_ary = zeros(nn, 1);
xi_ary(nn) = bb*exp(aalpha/sqrt(0.5 + H)*sqrt(N));
xi_ary(1:nn-1) = aa*(xi_ary(nn)/aa).^((1:(nn-1))/nn); 

cH = 1/(gamma(H+0.5)*gamma(0.5-H));

x_mt = zeros(nn, mm);
tilw_mt = zeros(nn, mm);

[x_mt(1, :), tilw_mt(1, :)] =   gauss_quadrature(mm, xi_ary(1), [], H, cH, 'Fraction');

for i = 2:nn
    [x_mt(i, :), tilw_mt(i, :)] =   gauss_quadrature(mm, xi_ary(i), xi_ary(i-1), [], [], 'Legendre');
end

w_mt = tilw_mt;
w_mt(2:end, :) = cH*x_mt(2:end, :).^(-H-0.5).*tilw_mt(2:end, :);

KNt = zeros(length(t_ary), 1);
for i = 1:length(t_ary)
    KNt(i) = sum(w_mt.*exp(-x_mt*t_ary(i)), 'all');
end
end