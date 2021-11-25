function [K] = compute_Kn0_3D(s, h, x1 , x2 , x3 , pts , wts)
% compute K_i for the cells not centered at the origin
[X1,~,~] = meshgrid(x1,x2,x3);
K = 0.*X1; % Initialize K matrix
for i3 = 1:length(x3) % loop over cells in x3 direction
    for i2 = 1:length(x2) % loop over cells in x2 direction
        for i1 = 1:length(x1) % loop over cells in x1 direction
            Ki = 0; % Initialize quadrature contributions for cell i1,i2,i3
            for q3 = 1:length(pts) % loop over quadrature points in x3 direction
                xq3 = x3(i3) + h*pts(q3)/2; % quadrature point x3 location (physical domain)
                for q2 = 1:length(pts) % loop over quadrature points in x2 direction
                    xq2 = x2(i2) + h*pts(q2)/2; % quadrature point x2 location (physical domain)
                    for q1 = 1:length(pts) % loop over quadrature points in x1 direction
                        xq1 = x1(i1) + h*pts(q1)/2; % quadrature point x1 location (physical domain)
                        Ki = Ki + pi^(-3/2)*2^(-2*s(xq1,xq2,xq3)) * (gamma(3/2 - s(xq1,xq2,xq3))/gamma(s(xq1,xq2,xq3)))...
                            *(1/(sqrt(xq1^2 + xq2^2 + xq3^2)^(3-2*s(xq1,xq2,xq3))))...
                            *(h/2)^3*(wts(q1)*wts(q2)*wts(q3))./(h^3); % quadrature point contribution 
                    end
                end
            end
            K(i1,i2,i3) = Ki;
        end 
    end
end
K(1,1,1) = 0;
end