function phi = char_var(params, tau, u, t_num, T, model, x_ary, w_ary)
% calculating the characteristic functions of integrated variance component for different models.
    if strcmp(model, 'rHeston_var')
        V1 = params(1);
        H = params(4);
        k1 = params(5);
        theta1 = params(6);
        sigma1 = params(7);
        alpha = H + 0.5;  
        u_num = size(u, 2);

        NN = length(x_ary);
        dt = T/t_num;   
        t_ary  = (0:(t_num-1)) / t_num * T;
        T_ary = T:dt:(T+tau);
        bb = -k1;
        g0_ary = V1 + k1*theta1/gamma(alpha+1)*t_ary.^alpha;
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
        integ1 = g0T_ary*(E1*bb + 1)/tau*dt * 1i*u;
        
        F_mt = zeros(t_num, u_num)*nan;
        initial_ary = 1;
        N_tmp = 10;
        while initial_ary(end)>1e-6
            n_ary = 0:N_tmp;
            initial_ary = bb.^n_ary.*tau.^(alpha*(n_ary+1)-1)./gamma(alpha*(n_ary+1)+1);
            N_tmp = N_tmp * 10;
            if N_tmp >= 1e5
               fprintf('<<Errors in calculating rHeston_var');
               break;
            end
        end
        for u_index = 1:u_num
            phi_ary = zeros(t_num, NN);
            phi_ary(1, :) = 1i*u(u_index)*sum(initial_ary);
            for t_index = 2:t_num+1
                sum_tmp = w_ary' * phi_ary(t_index-1, :).';
                F_mt(t_index-1, u_index) = -k1*sum_tmp + 1/2*sigma1^2*sum_tmp^2;
                phi_ary(t_index, :) = (phi_ary(t_index-1, :) + dt*F_mt(t_index-1))./(1+x_ary'*dt);
            end
        end
        integ2 = flip(g0_ary)*F_mt*dt;
        phi = exp(integ1 + integ2);


   elseif  strcmp(model,'Heston_var')
        V = params(1);
        k  = params(2);
        theta    = params(3);
        sigma     = params(4);       
        C = (1+1i*u*sigma^2/(2*k)*(exp(-k*T - 1))).^(-2*k*theta/sigma^2);
        D = exp(2*k*1i*u*V./(sigma^2*u*1i + (2*k - sigma^2*u*1i)*exp(k*T)));
        phi = C.*D;
    end

end

