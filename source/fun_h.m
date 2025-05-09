function dhdt = fun_h(y, h, k, sigma)

dhdt = zeros(2,1);

dhdt(1) = h(2);
dhdt(2) = 2/(sigma^2*y^2)*((k*y+1)*h(2) - 1);
end