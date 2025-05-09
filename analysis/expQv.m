function [exQV] = expQv(model, params, du)
% calculate the expected integrated variance under the Q measure.
tau = 1/12;

if strcmp(model, '42rswithmu')
    
    V1 = params(1); 
    H = params(4);
    k1 = params(5);
    theta1 = params(6);
    sigma1 = params(7);
    alpha = H + 0.5;  
    
    dt = tau/1e2;
    T = 0;
    T_ary = T:dt:(T+tau);
    bb = -k1;

    g0T_ary = V1 + k1*theta1/gamma(alpha+1)*T_ary.^alpha;
    E1 = zeros(length(T_ary), 1);
    for t_index = 1:length(T_ary)
        t_tmp = T_ary(t_index);
        ttmp = 1;
        N_tmp = 10;
        while abs(ttmp(end))>1e-6
            n_ary = 0:N_tmp;
            ttmp = bb.^n_ary./gamma(alpha*(n_ary+1)+1).*(tau - t_tmp + T).^(alpha*(n_ary + 1));
            N_tmp = N_tmp * 10;
            if N_tmp >= 1e5
                fprintf('<<Errors in calculating rHeston_var');
                break;
            end
        end
        E1(t_index) = sum(ttmp);    
    end
    expbarV1 = g0T_ary*(E1*bb + 1)/tau*dt;


    V2 = params(2);
    mu = params(3);
    k2 = params(9);
    theta2 = params(10);
    sigma2 = params(11);

    deltal = du;
    tilalpha = -(0.5 + k2/sigma2^2) + sqrt((0.5 + k2/sigma2^2).^2 + 2*deltal/sigma2^2);
    tilgam = 2*(tilalpha + 1 + k2/sigma2^2);
    tily = V2/(k2*theta2)*(exp(k2*theta2*tau)-1);
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
    expbarV2 = h_ary/tau;
    exQV = mu*expbarV1 + (1-mu)*expbarV2;


elseif strcmp(model, '42withmu')
    V1 = params(1); 
    k1 = params(4);
    theta1 = params(5);
    sigma1 = params(6);
    
    aa1 = (1-exp(-k1*tau))/(k1*tau);
    bb1 = theta1*(1-aa1);
    expbarV1 = aa1*V1 + bb1;

    V2 = params(2);
    mu = params(3);
    k2 = params(8);
    theta2 = params(9);
    sigma2 = params(10);

    deltal = du;
    tilalpha = -(0.5 + k2/sigma2^2) + sqrt((0.5 + k2/sigma2^2).^2 + 2*deltal/sigma2^2);
    tilgam = 2*(tilalpha + 1 + k2/sigma2^2);
    tily = V2/(k2*theta2)*(exp(k2*theta2*tau)-1);
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
    expbarV2 = h_ary/tau;

    exQV = mu*expbarV1 + (1-mu)*expbarV2;

elseif strcmp(model, 'H2Fwithmu')
    mu = params(3);
    V1 = params(1); 
    V2 = params(2); 
    k1 = params(4);
    theta1 = params(5);
    sigma1 = params(6);
    
    k2 = params(8);
    theta2 = params(9);
    sigma2 = params(10);

    aa1 = (1-exp(-k1*tau))/(k1*tau);
    bb1 = theta1*(1-aa1);
    expbarV1 = aa1*V1 + bb1;

    aa2 = (1-exp(-k2*tau))/(k2*tau);
    bb2 = theta2*(1-aa2);

    expbarV2 = aa2*V2 + bb2;
    exQV = mu*expbarV1 + (1-mu)*expbarV2;    

end
end

