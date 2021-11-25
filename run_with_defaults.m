% Run full implementation with default settings
s = @(x1,x2,x3) 0.1.*(6 + 9.*sqrt(x1.^2+x2.^2+x3.^2))/(1+ sqrt(x1.^2+x2.^2+x3.^2)); % standard (3D) gravitational field
%s = @(r) 0.1.*(6 + 9.*r)/(1+r); % standard (radial) gravitational field
h = 0.2; % cell size
R = 1; % mass distribution (source) radius
D = 5; % domain radius (K domain will have double this radius)
K_opt = "3D"; % option to find FK in 3D or radially ("3D" or "radial")
plot_K = true; % (boolean) option to plot K (modified gravitational field) on
plot_I = true; % (boolean) option to plot I (gravitational potential)
plot_v_orb = true; % (boolean) option to plot v_orb (orbital velocity)
compare_to_exact = false; % boolean (option to compare to the exact solution (only for constant s)
[x1, x1K, K, f, I_3D, err_I] = full_implementation(s, h, R, D, K_opt, plot_K, plot_I, plot_v_orb, compare_to_exact);