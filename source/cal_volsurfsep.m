function [error] = cal_volsurfsep(params, spx_tb, vix_tb, model, method, opt)
% calculating SPX and VIX option prices seperately.
% SPX option pricing
if isempty(spx_tb)
    diff_ary = [];
else
    M = size(spx_tb,1);
    expiries = unique(spx_tb.t,'stable');
    L = size(expiries,1);
    Lcont = 10;
    c1 = 0.1*expiries;
    c2 = 0.04*expiries;
    Xh = c1+sqrt(c2)*Lcont;
    Xl = c1-sqrt(c2)*Lcont;
    Xc = (Xh-Xl)/2;
    Xm = (Xh+Xl)/2;
    
    counter = 1:L;
    IR = spx_tb.interest_rate(1);
    DY = spx_tb.dividend_rate(1);
    power = 10;
    grid_num = 100;
    CP = 1;

    diff_ary = zeros(M, 1)*nan;
    count = 1;
    for i = 1:length(counter)
        expiry = expiries(i);
        pick = spx_tb.t == expiry;
        tmp = spx_tb(pick,:);
        T    = unique(tmp.t);
        K    = tmp.K;
        m    = length(K);
        S0 = tmp.S(1);
        vg = tmp.vega;
        option_price = tmp.mid_price;
        if strcmp(method, 'COS')
            smile_tmp = COS_pricing(S0,T,K,IR,DY,model,params,grid_num,Xc(i),Xm(i),2^power,CP,1);
        end
        if all(isnan(smile_tmp))
            smile_tmp(isnan(smile_tmp)) = 0;
        else
            smile_tmp(isnan(smile_tmp)) = mean(smile_tmp, 'omitnan') ;
        end
        diff_ary(count: count+m-1, 1) = (smile_tmp -option_price)./vg;
        count = count + m;
    end
end

% VIX option pricing
if isempty(vix_tb)
    diff2_ary = [];
else
    IR = vix_tb.interest_rate(1);
    MM = size(vix_tb,1);
    expiries2 = unique(vix_tb.t,'stable');
    expiries2 = sort(expiries2);
    L2 = size(expiries2,1);
    tau = 30/365;
    count = 1;
    diff2_ary = zeros(MM, 1)*nan;
    upV = 1.5;
    if year(vix_tb.date(1)) < 2017 || year(vix_tb.date(1)) > 2018
        upV = upV + 1;  
        % ``[0, upV]" is the range of the integration of density functions
        %  and is adjusted for all the models to ensure the density is not seriously truncated. 
    end

    if strcmp(model, '42rswithmu') && params(4) < 0.4999
        % calculating function h2 in Prop 3.3.
        tgrid_num = 1e3;
        N = 2^10;
        charrough = cal_char42rs(params, tau, N, 0, upV, tgrid_num, expiries2);
        x_num = 1e3;
        x_density = linspace(0, upV, x_num);
        mu = params(3);
        a = sqrt(mu);
        V2 = params(2);
        b = sqrt(1-mu);
        k2 = params(9);
        theta2 = params(10);
        sigma2 = params(11);
        y_num = 1e3;
        y_density = linspace(1e-4, upV, y_num);
        deltal = 1e-8;
        tilalpha = -(0.5 + k2/sigma2^2) + sqrt((0.5 + k2/sigma2^2).^2 + 2*deltal/sigma2^2);
        tilgam = 2*(tilalpha + 1 + k2/sigma2^2);
        tily = y_density/(k2*theta2)*(exp(k2*theta2*tau)-1);
        prod_M11 = 1;
        NN = 10;
        tilz = 2./(sigma2^2*tily);
        while max(abs(prod_M11(end, :)))>1e-6
            elem_M11 = (tilalpha + (0:1:(NN-1))')./(tilgam + (0:1:(NN-1))')*(-1*tilz)./(1:1:NN)';
            prod_M11 = cumprod(elem_M11, 1);
            M11 = 1 + sum(prod_M11);
            NN = NN*10;
            if NN >= 1e5
                fprintf('<<Errors in calculating M11');
                break;
            end
        end
        e_inte = gamma(tilgam - tilalpha)/gamma(tilgam)*tilz.^tilalpha.*M11;
        h_ary = -(e_inte-1)/deltal;
        h_ary(h_ary<0) = 0;
   end
    
    for i = 1:L2
        T = expiries2(i);
        pick = vix_tb.t == T;
        tmp = vix_tb(pick,:);
        K    = tmp.K;
        m    = length(K);
        option_price = tmp.mid_price;
        vg = tmp.vega;
        if strcmp(model, '42rswithmu') && params(4) < 0.4999
              f1_ary = COS_roughdensity(x_density, 0, upV, N, charrough(i, :));
              f2_ary = density32(y_density, T, k2, sigma2, theta2, V2);
              dx = x_density(2) - x_density(1);
              dy = y_density(2) - y_density(1);
              vix_price = zeros(length(K), 1);
              for j = 1:length(K)
                  mt = max(sqrt(a^2*x_density + b^2/tau*h_ary')*100- K(j), 0);% !!! 
                  vix_price(j) = f2_ary*mt*f1_ary*dy*dx;
              end
              vix_price = vix_price*exp(-IR*T);
        elseif strcmp(model, '42rswithmu') && params(4) >= 0.4999
            params_tmp = [params(1:3); params(5:end)];
            vix_price = cal_vixprice(params_tmp, tau, K, T, IR, '42withmu');
        else 
            vix_price = cal_vixprice(params, tau, K, T, IR, model);
        end
        
        if all(isnan(vix_price))
            vix_price(isnan(vix_price)) = 0;
        else
            vix_price(isnan(vix_price)) = mean(vix_price, 'omitnan') ;
        end

        diff2_ary(count: count+m-1, 1) = (vix_price - option_price)./vg;
        count = count + m;
    end

end

error_all = [diff_ary; diff2_ary];
if strcmp(opt,'ip')
    error = sum(error_all.^2);
elseif strcmp(opt,'lm')
    error = error_all;
end
end
