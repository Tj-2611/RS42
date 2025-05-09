function f_ary = COS_density(x, T, tau, model, params, t_num, a, b, N, x_ary, w_ary)
% calculating density function using cos method (Prop A.1.)
    k   = 0:N-1;
    u   = k*pi/(b-a);

    phi = char_var(params, tau, u, t_num, T, model, x_ary, w_ary);
    phi(1) = 1;
    m   = length(x);
    eto = zeros(m,N);
    for count = 1:m
        eto(count,:) = exp(1i*pi*k*(-a)/(b-a));
    end
    
    wgt = ones(1,N);
    wgt(1) = 0.5;
    
    f_ary = zeros(m,1)*nan;
    for count = 1:m
        psi = cos(k*pi*(x(count)-a)/(b-a));
        f_ary(count) = 2/(b-a)*real(wgt.*phi.*eto(count,:))*psi';
    end
    
    
end