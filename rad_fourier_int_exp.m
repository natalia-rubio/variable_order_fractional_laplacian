function [FK, dK_vec] = rad_fourier_int_exp(k_val, lambda, s, max_iter, domain_size)
    K_s = @(r) gamma(1.5 - s(r))./(4.^s(r) .* pi.^(3/2) .* gamma(s(r))...
        .* r.^(2.*(1.5 - s(r)))); % specify K_s
    T = 2*pi / k_val ; % period of integration
    i_T = 0; % integration iteration counter
    cvg_flag = false; % convergence flag
    FK = 0; % initialize FK
    dK_vec = zeros(max_iter, 1); % dK integral increment vector

    if k_val == 0 % k = 0 case
       while cvg_flag == false % loop over integration increments
            T = 1; % set integration period
            fun = @(r) exp(-lambda.*r) .* r .* K_s(r) .* (r); % set integrand
            dK = (4*pi  * integral(fun, T*i_T, T*(i_T+1))); % find integration increment
            FK = FK + dK; % add increment contribution

            cvg_flag = ((abs(dK/FK) + abs(dK) < 1e-7) | (1+i_T >= max_iter)); % check convergence flag
            i_T = i_T + 1; % increase integration increment counter
            dK_vec(i_T) = dK; % store integration increment in vector
           
            if i_T == max_iter % in case of non-convergence, display relevant values
                disp("failure to converge for k_val = ")
                k_val
                FK
                dK
            end
            if T*(i_T+1) > domain_size % restrict integration to the domain
                cvg_flag=true;
            end
       end
    else
        while cvg_flag == false
            fun = @(r) exp(-lambda.*r) .* r .* K_s(r) .* sin(2*pi*k_val.*r); % set integrand
            dK = (4*pi/(2*pi*k_val)) * integral(fun, T*i_T, T*(i_T+1)); % find integration increment 
            FK = FK + dK; % add increment contribution

            cvg_flag = ((abs(dK/FK) + abs(dK) < 1e-7) | (1+i_T >= max_iter)); % check convergence flag
            i_T = i_T + 1; % increase integration increment counter
            dK_vec(i_T) = dK; % store integration increment in vector
           
            if i_T == max_iter % in case of non-convergence, display relevant values
                disp("failure to converge for k_val = ")
                k_val
                FK
                dK
            end
            if T*(i_T+1) > 100*domain_size % restrict integration to the domain
                cvg_flag=true;
            end
        end
    end
    dK_vec = dK_vec(1:i_T); % trim dK_vector
end