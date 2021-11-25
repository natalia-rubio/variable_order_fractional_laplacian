function [FK_radial_3D, K_radial_3D] = map_radial_to_3D(FK_radial, x1K, x2K, x3K, K_s)
[X1,X2,X3] = meshgrid(x1K,x2K,x3K); % generate 3D grid

spline = csaps(linspace(min(x1K),max(x1K), length(FK_radial)), FK_radial,1); % generate spline
spline_x = fnxtr(spline,2); % generate extrapolation function
FK_radial_3D = fnval(spline_x,sqrt(X1.^2 + X2.^2 + X3.^2));

K_radial_3D = K_s(sqrt(X1.^2 + X2.^2 + X3.^2)); % evaluate K at 3D grid locations
K_radial_3D(1,1,1) = 0; % set origin element to zero
end

