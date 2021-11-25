% Computation of K (modified gravitational field), f (mass distribution
% source), and I (gravitational potential).  Option to plot K, I, and orbital
% velocities.  Option to compare to exact solution (only for constant s)


% INPUTS:
% s --> function handle (function of x1,x2,x3 Cartesian coordinates if
% K_opt = "3D" or function of r radial coordinate if K_opt = "radial")
% h --> element size
% R --> source radius
% D --> domain radius
% plot_K --> (boolean) option to plot K (modified gravitational field) 
% plot_I --> (boolean) option to plot I (gravitational potential)
% plot_v_orb --> (boolean) option to plot v_orb (orbital velocity)
% compare_to_exact --> (boolean) option to compare to the exact solution (only for
% constant s)

% OUTPUTS:
% x1 --> spatial discretization for f, I (0-centered)
% x1K --> spatial discretization for K (0-centered)
% K --> computed K values
% f --> computed f (mass distribution) values
% I --> computed I (potential) values
% err_I --> relative error compared to explicit solution (only for constant
% s)

function [x1, x1K, K, f, I, err_I] = full_implementation(s, h, R, D, ...
    K_opt, plot_K, plot_I, plot_v_orb, compare_to_exact)
%% PLOTTING STYLE:
set(0,'defaulttextInterpreter','latex')
set(0,'defaultlegendInterpreter','latex')
set(0, 'defaultAxesFontName', 'Times New Roman')
set(0, 'DefaultLineLineWidth', 2);
        
%% SET UP MESH
num_el = ceil(D/h); % 1/2 number of elements on side length of domain
x1 = linspace(0 , h*(num_el-1), num_el); % x1 locations of cell centers
x2 = linspace(0 , h*(num_el-1), num_el); % x2 locations of cell centers
x3 = linspace(0 , h*(num_el-1), num_el); % x3 locations of cell centers

num_elK = 2*(num_el-1)+1; % double the size of the K domain for convolution
x1K = linspace(0 , h*(num_elK-1), num_elK); % x1 locations of K cell centers
x2K = linspace(0 , h*(num_elK-1), num_elK); % x2 locations of K cell centers
x3K = linspace(0 , h*(num_elK-1), num_elK); % x3 locations of K cell centers

domain_size = max(x1K)+h/2; % K domain radius
K_ramp_filter = oct2full(K_ramp(x1K,x2K,x3K)); % filter to eliminate destabilizing boundary effects

%% COMPUTE f_i:
num_quad = 30; % number of quadrature points for calculation of f_i
f = compute_f(h, R, x1 , x2 , x3, num_quad); % compute the f_i
f = oct2full(f); % convert octant to full domain
f = f.*(-3*R^-3); % normalize by source volume

%% FIND K & EXECUTE CONVOLUTION
if K_opt == "3D"
    
    % COMPUTE K_i in 3D:
    num_quad_K = 3; % number of quadrature points for calculation of K_i
    [pts, wts] = lgwt(num_quad_K,-1,1); % generate quadrature points for calculation of K_i
    K_3D = compute_Kn0_3D(s, h, x1K, x2K, x3K, pts, wts); % compute the K_i for cell not containing the origin
    K0_3D = compute_K01_3D(s, h, pts , wts) + compute_K02_3D(s, h, pts, wts); % compute K_0 for the cell centered at the origin
    K_3D(1,1,1) = K0_3D/h^3; % replace the origin centered cell
    K_3D = oct2full(K_3D); % convert octant to full domain
    K = K_3D; % set K
    % COMPUTE I_i:
    [I] = convolve_3D_K(K_3D.*K_ramp_filter,f); % convolve K_i and f_i
    
    if plot_K == true
        figure(1)
        dk = 1/((2*num_elK-1)*h); % Fourier space resolution
        k = dk*(0:((2*num_elK-1)/2)); % Fourier space frequencies
        FK = fftn(ifftshift(K_3D.*K_ramp_filter))./((2*num_elK-1)^3); % compute Fourier transform of filtered K_i
        FK_line = FK(1:num_elK,1,1); % extract radial line
        subplot(1,2,1)
        semilogy(k,abs(FK_line))
        title("$\widehat{K}_{s(\cdot)} (\textbf{x}) (\textbf{k})$ (3D)")
        xlabel("$\textbf{k}$"); ylabel("$\widehat{K}_{s(\cdot)}$")
        set(gca,'FontSize',14); hold off

        subplot(1,2,2)
        semilogy(x1K,-1.* K_3D(num_elK:end,num_elK,num_elK)); hold on
        title("$\Phi = - {K}_{s(\cdot)}(r=|\textbf{x}|)$")
        xlabel("$|\textbf{x}|$"); ylabel("$- {K}_{s(\cdot)}$")
        set(gca,'FontSize',14); hold off
    end
        % PLOT I:
    if plot_I == true
        figure(2)
        surf([-1.*flip(x1),x1(2:end)],[-1.*flip(x2),x2(2:end)],I(:,:,num_el)) % plot I on a center plane
        title(sprintf("Gravitational Potential (3D)"))
        xlabel("$x_1$"); ylabel("$x_2$");zlabel("potential")
        shading interp
        set(gca,'FontSize',14); hold off
    end

    % PLOT ORBITAL VELOCITIES:
    if plot_v_orb == true
        figure(3)
        pot_vals = I(num_el:end,num_el,num_el); % potential values on a radial line
        v_orb = zeros(length(pot_vals)-1,1); % initialize orbital velocity vector
        v_orb_locs = x1(1:end-1) + h/2; % r locations
        for i = 2:length(pot_vals)
            v_orb(i-1) = (pot_vals(i)-pot_vals(i-1))*v_orb_locs(i-1)/h; % finite difference for orbital velocities
        end
        plot(v_orb_locs, v_orb) % plot orbital velocities along a radial line
        title(sprintf("Orbital Velocity (3D)"))
        xlabel("$r$"); ylabel("$\Omega = r \frac{dI}{dr}$")
        set(gca,'FontSize',14); hold off
    end

    % COMPARE TO EXACT SOLUTION:
    if compare_to_exact == true
        [I_exp] = explicit_sol(x1,x2,x3, R, s(0,0,0)); % compute explicit solution
        I_exp = oct2full(I_exp); % convert octant to full domain
        err_I = abs(I-I_exp)./I_exp; % relative error between numerical and explicit Is
    elseif compare_to_exact == false
        err_I = NaN;
    end 
elseif K_opt == "radial"

    % COMPUTE K_i RADIALLY:
    K_s = @(r) gamma(1.5 - s(r))./(4.^s(r) .* pi.^(3/2) .* gamma(s(r))...
        .* r.^(2.*(1.5 - s(r)))); % specify K_s
    K = K_s(x1K); % set K
    [FK_radial, ~] = compute_FK_radial(s, x1K, h, domain_size); % Compute FK radially
    [FK_radial_3D, ~] = map_radial_to_3D(FK_radial, x1K, x2K, x3K, K_s); % Map radial FK into 3D

    % COMPUTE I_i:
    [I] = convolve_FK_radial(FK_radial_3D,f,domain_size); % convolve radial FK_i and f_i
    %I(num_el,num_el,num_el) = NaN; % set origin cell to NaN
    if plot_K == true
        FK_radial_line = FK_radial_3D(:,1,1)/((2*domain_size)^3);
        figure(4); subplot(1,2,1); semilogy(x1K, FK_radial_line); hold on; subplot(1,2,2)
        semilogy(x1K, -K_s(x1K)); hold on; subplot(1,2,1)
        title("$\widehat{K}_{s(\cdot)} (\textbf{x}) (\textbf{k})$ (radial)")
        xlabel("$\textbf{k}$"); ylabel("$\widehat{K}_{s(\cdot)}, \; \lambda=0.1$")

        set(gca,'FontSize',14); hold on; subplot(1,2,2)
        title("$\Phi = - {K}_{s(\cdot)}(r=|\textbf{x}|)$")
        xlabel("$|\textbf{x}|$"); ylabel("$- {K}_{s(\cdot)}$")
        set(gca,'FontSize',14); hold off
    end
            % PLOT I:
    if plot_I == true
        figure(5)
        surf([-1.*flip(x1),x1(2:end)],[-1.*flip(x2),x2(2:end)],I(:,:,num_el)) % plot I on a center plane
        title(sprintf("Gravitational Potential (radial)"))
        xlabel("$x_1$"); ylabel("$x_2$");zlabel("potential")
        shading interp
        set(gca,'FontSize',14); hold off
    end
    % PLOT ORBITAL VELOCITIES:
    if plot_v_orb == true
        figure(6)
        pot_vals = I(num_el:end,num_el,num_el); % potential values on a radial line
        v_orb = zeros(length(pot_vals)-1,1); % initialize orbital velocity vector
        v_orb_locs = x1(1:end-1) + h/2; % r locations
        for i = 2:length(pot_vals)
            v_orb(i-1) = (pot_vals(i)-pot_vals(i-1))*v_orb_locs(i-1)/h; % finite difference for orbital velocities
        end
        plot(v_orb_locs, v_orb) % plot orbital velocities along a radial line
        hold on
        title(sprintf("Orbital Velocity (radial)"))
        xlabel("$r$"); ylabel("$\Omega = r \frac{dI}{dr}$")
        set(gca,'FontSize',14); hold off
    end
    % COMPARE TO EXACT SOLUTION:
    if compare_to_exact == true
        [I_exp] = explicit_sol(x1,x2,x3, R, s(0,0,0));
        I_exp = oct2full(I_exp); % convert octant to full domain
        err_I = abs(I-I_exp)./I_exp; % relative error between numerical and explicit Is
    elseif compare_to_exact == false
        err_I = NaN;
    end
else
        print("K_opt must be '3D' or 'radial'.")
end
end
