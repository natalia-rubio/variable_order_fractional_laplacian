function [K02] = compute_K02_3D(s, h, pts , wts)
% compute k_2 contribution for cell centered at origin
K02 = 0; % Initialize quadrature contributions for cell
for q = 1:length(pts) % loop over quadrature points in radial direction
    rhoq = h/4 + h*pts(q)/4; % quadrature point rho location (physical domain)
    K02 = K02 + (4*pi)*pi^(-3/2)*2^(-2*s(rhoq,0,0)) * gamma(3/2 - s(rhoq,0,0))/gamma(s(rhoq,0,0))...
    *(1-ramp_func(rhoq*2/h))*(1/(rhoq^(1-2*s(rhoq,0,0))))*(h/4)*(wts(q)); % quadrature point contribution
end
end