function [vix_price] = cal_vixprice(params, tau, K, T, IR, model)
% calculating vix option prices for different models
    if strcmp(model, '42rswithmu')
        mu = params(3);
        a = sqrt(mu);
        H = params(4);
        V2 = params(2);
        b = sqrt(1-mu);
        tgrid_num = 1e3;
        x_num = 1e3;
        upV = T + 1.5;
        x_density = linspace(0, upV, x_num);
        N = 2^10;
        dt = T/tgrid_num;
        [xn_ary, wn_ary] = markovianappr(H, T, dt, tgrid_num, 20);
        f1_ary = COS_density(x_density, T, tau, 'rHeston_var', params, tgrid_num, 0, upV, N, xn_ary, wn_ary);
        k2 = params(9);
        theta2 = params(10);
        sigma2 = params(11);
        y_num = 5e3;
        y_density = linspace(1e-4, upV, y_num);
        f2_ary = density32(y_density, T, k2, sigma2, theta2, V2);
        deltal = 1e-8;
        tilalpha = -(0.5 + k2/sigma2^2) + sqrt((0.5 + k2/sigma2^2).^2 + 2*deltal/sigma2^2);
        tilgam = 2*(tilalpha + 1 + k2/sigma2^2);
        tily = y_density/(k2*theta2)*(exp(k2*theta2*tau)-1);
        prod_M11 = 1;
        N = 10;
        tilz = 2./(sigma2^2*tily);
        while max(abs(prod_M11(end, :)))>1e-6
            elem_M11 = (tilalpha + (0:1:(N-1))')./(tilgam + (0:1:(N-1))')*(-1*tilz)./(1:1:N)';
            prod_M11 = cumprod(elem_M11, 1);
            M11 = 1 + sum(prod_M11);
            N = N*10;
            if N >= 1e5
                fprintf('<<Errors in calculating M11')
                break;
            end
        end
        e_inte = gamma(tilgam - tilalpha)/gamma(tilgam)*tilz.^tilalpha.*M11;
        h_ary = -(e_inte-1)/deltal;
        h_ary(h_ary<0) = 0;

        dx = x_density(2) - x_density(1);
        dy = y_density(2) - y_density(1);
        vix_price = zeros(length(K), 1);
        for j = 1:length(K)
            mt = max(sqrt(a^2*x_density + b^2/tau*h_ary')*100- K(j), 0);% !!! 
            vix_price(j) = f2_ary*mt*f1_ary*dy*dx;
        end
        vix_price = vix_price*exp(-IR*T);
    elseif strcmp(model, '42withmu')
        V1 = params(1);
        mu = params(3);
        a = sqrt(mu);
        b = sqrt(1-mu);
        params1 = [V1; params(4:7)];
        x_num = 1e3;
        upV = T + 1.5;
        x_density = linspace(0, upV, x_num);
        N = 2^10;       
        f1_ary = density_Heston(x_density, T, params(4), params(6), params(5), params(1));
        V2 = params(2);       
        k2 = params(8);
        theta2 = params(9);
        sigma2 = params(10);
        y_num = 1e3;
        y_density = linspace(1e-4, upV, y_num);
        f2_ary = density32(y_density, T, k2, sigma2, theta2, V2);
        deltal = 1e-8;
        tilalpha = -(0.5 + k2/sigma2^2) + sqrt((0.5 + k2/sigma2^2).^2 + 2*deltal/sigma2^2);
        tilgam = 2*(tilalpha + 1 + k2/sigma2^2);
        tily = y_density/(k2*theta2)*(exp(k2*theta2*tau)-1);
        prod_M11 = 1;
        N = 10;
        tilz = 2./(sigma2^2*tily);
        while max(abs(prod_M11(end, :)))>1e-6
            elem_M11 = (tilalpha + (0:1:(N-1))')./(tilgam + (0:1:(N-1))')*(-1*tilz)./(1:1:N)';
            prod_M11 = cumprod(elem_M11, 1);
            M11 = 1 + sum(prod_M11);
            N = N*10;
            if N >= 1e5
                fprintf('<<Errors in calculating M11')
                break;
            end
        end
        e_inte = gamma(tilgam - tilalpha)/gamma(tilgam)*tilz.^tilalpha.*M11;
        h_ary = -(e_inte-1)/deltal;
        h_ary(h_ary<0) = 0;
        
        dx = x_density(2) - x_density(1);
        dy = y_density(2) - y_density(1);
             
        k1 = params1(2);
        theta1 = params1(3);
        aa1 = (1-exp(-k1*tau))/(k1*tau);
        bb1 = theta1*(1-aa1);
        vix_price = zeros(length(K), 1);
        for j = 1:length(K)
            mt = max(sqrt(a^2*(aa1*x_density + bb1) + b^2/tau*h_ary')*100- K(j), 0);
            vix_price(j) = f2_ary*mt*f1_ary*dy*dx;
        end
        vix_price = vix_price*exp(-IR*T);        

    elseif strcmp(model, 'H2Fwithmu')
        V1 = params(1);
        V2 = params(2);
        mu = params(3);
        a = sqrt(mu);
        b = sqrt(1-mu);
        params1 = [V1; params(4:7)];
        params2 = [V2; params(8:11)];
        x_num = 1e3;
        upV = T + 1.5;
        x_density = linspace(0, upV, x_num);
        N = 2^10;
        f1_ary = density_Heston(x_density, T, params(4), params(6), params(5), params(1));
        f2_ary = density_Heston(x_density, T, params(8), params(10), params(9), params(2));
        dx = x_density(2) - x_density(1);
        vix_price = zeros(length(K), 1);
        k1 = params1(2);
        k2 = params2(2);
        theta1 = params1(3);
        theta2 = params2(3);
        aa1 = (1-exp(-k1*tau))/(k1*tau);
        bb1 = theta1*(1-aa1);
        aa2 = (1-exp(-k2*tau))/(k2*tau);
        bb2 = theta2*(1-aa2);
        for j = 1:length(K)
            mt = max(sqrt(a^2*(aa1*x_density + bb1) + b^2*(aa2*x_density' + bb2))*100- K(j), 0);
            vix_price(j) = f2_ary'*mt*f1_ary*dx*dx;
        end
        vix_price = vix_price*exp(-IR*T);    
    end
end
