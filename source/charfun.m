function phi = charfun(model,params,n,u,t)
% calculating the characteristic functions of the log price for different models.
    if strcmp(model,'H2Fwithmu')
       mu = params(3);
        a = sqrt(mu);
        b = sqrt(1-mu);
       V1 = params(1)*a^2; 
       V2 = params(2)*b^2; 
        lambda = params(4);
        v_bar  = params(5)*a^2;
        eta    = params(6)*a;
        rho    = params(7);
        beta  = lambda - rho*eta*1i*u;
        alpha = -0.5*u.^2 - 0.5*1i*u;
        zeta  = 0.5*eta^2;
        d = sqrt(beta.^2-4*alpha*zeta);
        rm = (beta-d)/eta^2;
        rp = (beta+d)/eta^2;
        g = rm./rp;

        C = rm*t - 2/eta^2*log((1-g.*exp(-d*t))./(1-g));
        D = rm.*(1-exp(-d*t))./(1-g.*exp(-d*t));

        phi1 = exp(v_bar*lambda*C + V1*D);
        
        lambda2 = params(8);
        v_bar2  = params(9)*b^2;
        eta2    = params(10)*b;
        rho2    = params(11);
        beta2  = lambda2 - rho2*eta2*1i*u;
        alpha2 = -0.5*u.^2 - 0.5*1i*u;
        zeta2  = 0.5*eta2^2;
        d2 = sqrt(beta2.^2-4*alpha2*zeta2);
        rm2 = (beta2-d2)/eta2^2;
        rp2 = (beta2+d2)/eta2^2;
        g2 = rm2./rp2;

        C2 = rm2*t - 2/eta2^2*log((1-g2.*exp(-d2*t))./(1-g2));
        D2 = rm2.*(1-exp(-d2*t))./(1-g2.*exp(-d2*t));

        phi2 = exp(v_bar2*lambda2*C2 + V2*D2);     
        phi = phi1.*phi2;
    elseif strcmp(model, '42withmu')
        mu = params(3);
        a = sqrt(mu);
        b = sqrt(1-mu);
        V1 = params(1)*a^2; 
        lambda = params(4);
        v_bar  = params(5)*a^2;
        eta    = params(6)*a;
        rho    = params(7);
        beta  = lambda - rho*eta*1i*u;
        alpha = -0.5*u.^2 - 0.5*1i*u;
        zeta  = 0.5*eta^2;
        d = sqrt(beta.^2-4*alpha*zeta);
        rm = (beta-d)/eta^2;
        rp = (beta+d)/eta^2;
        g = rm./rp;

        C = rm*t - 2/eta^2*log((1-g.*exp(-d*t))./(1-g));
        D = rm.*(1-exp(-d*t))./(1-g.*exp(-d*t));

        phi1 = exp(v_bar*lambda*C + V1*D);
        
        V2 = params(2);
        k2 = params(8);
        theta2 = params(9);
        sigma2 = params(10);
        rho2 = params(11);
        
        tmpp = -k2 + 1i*sigma2*rho2*u*b;
        tilalpha = -(0.5 - tmpp/sigma2^2) + sqrt((0.5 - tmpp/sigma2^2).^2 + b^2/sigma2^2*(1i*u + u.^2));
        tilgam = 2*(tilalpha + 1- tmpp/sigma2^2);
        
        prod_M11 = 1;
        N = 1e2;
        z =  2*k2*theta2/(sigma2^2*V2*(exp(k2*theta2*t)-1));
        while max(prod_M11(end, :))>1e-6
            elem_M11 = (tilalpha + (0:1:(N-1))')./(tilgam + (0:1:(N-1))')*(-1*z)./(1:1:N)';
            prod_M11 = cumprod(elem_M11, 1);
            M11 = 1 + sum(prod_M11);
            N = N*10;
            if N >= 1e5
                fprintf('<<Errors in calculating M11');
                break;
            end
        end
        char32 = double(gamma(sym(tilgam - tilalpha))./gamma(sym(tilgam)).*sym(z).^tilalpha).*M11;
        phi = phi1.*char32;
    elseif strcmp(model, '42rswithmu')
        mu = params(3);
        a = sqrt(mu);
        V1 = params(1)*a^2;
        H1 = params(4);
        k1 = params(5);
        theta1 = params(6)*a^2;
        sigma1 = params(7)*a;
        rho1 = params(8);        

        V2 = params(2);
        b = sqrt(1-mu);
        k2 = params(9);
        theta2 = params(10);
        sigma2 = params(11);
        rho2 = params(12);
        
        alpha1 = H1 + 0.5;
        ti  = (0:(n-1)) / n * t;
        y1   = ti.^alpha1;
        [rm,rp,p1,p2,p3,p4,q1,q2,q3,q4] = dh_Pade44_coeff(alpha1, k1, rho1, sigma1, u);
        m = size(u,2);
        Y10T = zeros(1,m);
        for j = 1:m
            [h_pade, dah] = dh_Pade44(y1, sigma1, rm(1,j),rp(1,j),p1(1,j),p2(1,j),p3(1,j),p4(1,j),q1(1,j),q2(1,j),q3(1,j),q4(1,j)); 
            Y10T(1, j) = k1*theta1*sum(h_pade)*t/n + V1*sum(dah)*t/n;
        end
        
        tmpp = -k2 + 1i*sigma2*rho2*u*b;
        tilalpha = -(0.5 - tmpp/sigma2^2) + sqrt((0.5 - tmpp/sigma2^2).^2 + b^2/sigma2^2*(1i*u + u.^2));
        tilgam = 2*(tilalpha + 1- tmpp/sigma2^2);
        
        prod_M11 = 1;
        N = 1e2;
        z =  2*k2*theta2/(sigma2^2*V2*(exp(k2*theta2*t)-1));
        while max(prod_M11(end, :))>1e-6
            elem_M11 = (tilalpha + (0:1:(N-1))')./(tilgam + (0:1:(N-1))')*(-1*z)./(1:1:N)';
            prod_M11 = cumprod(elem_M11, 1);
            M11 = 1 + sum(prod_M11);
            N = N*10;
            if N >= 1e5
                fprintf('<<Errors in calculating M11');
                break;
            end
        end
        char32 = double(gamma(sym(tilgam - tilalpha))./gamma(sym(tilgam)).*sym(z).^tilalpha).*M11;
        phi = exp(Y10T).*char32;

    end

end
