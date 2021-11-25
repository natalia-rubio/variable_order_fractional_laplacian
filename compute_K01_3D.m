function [K01] = compute_K01_3D(s, h, pts , wts)
% compute k_1 contribution for cell centered at origin
K01 = 0; % Initialize quadrature contributions for cell
for q3 = 1:length(pts) % loop over quadrature points in x3 direction
    xq3 = h/4 + h*pts(q3)/4; % quadrature point x3 location (physical domain)
    for q2 = 1:length(pts) % loop over quadrature points in x2 direction
        xq2 = h/4 + h*pts(q2)/4; % quadrature point x2 location (physical domain)
        for q1 = 1:length(pts) % loop over quadrature points in x1 direction
            xq1 = h/4 + h*pts(q1)/4; % quadrature point x1 location (physical domain)
            r = sqrt(xq1^2 + xq2^2 + xq3^2);
            K01 = K01 + 8*pi^(-3/2)*2^(-2*s(xq1,xq2,xq3)) * gamma(3/2 - s(xq1,xq2,xq3))/gamma(s(xq1,xq2,xq3))...
                *(ramp_func(2*r/h)/(sqrt(xq1^2 + xq2^2 + xq3^2)^(3-2*s(xq1,xq2,xq3))))*...
                (h/4)^3*(wts(q1)*wts(q2)*wts(q3)); % quadrature point contribution
        end
    end
end
end