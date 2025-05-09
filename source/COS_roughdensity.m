function f_ary = COS_roughdensity(x, a, b, N, phi)
% calculating the density function of integrated rough component (Prop A.1)
    k   = 0:N-1;
    u   = k*pi/(b-a);
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