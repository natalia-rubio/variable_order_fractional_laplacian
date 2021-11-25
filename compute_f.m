function [f] = compute_f(h, R, x1 , x2 , x3, num_quad)
% compute f over all cells
[pts, wts] = lgwt(num_quad,-1,1); % Set quadrature points and weights
[f,~,~] = meshgrid(x1,x2,x3);
f = 0.*f; % Initialize f

for i3 = 1:length(x3) % loop over elements in x3 direction
    for i2 = 1:length(x2) % loop over elements in x2 direction
        for i1 = 1:length(x1) % loop over elements in x1 direction
            fi = 0;
            r1 = sqrt((x1(i1)-h/2)^2 + (x2(i2)-h/2)^2 + (x3(i3)-h/2)^2); % radius of inner corner of cell
            r2 = sqrt((x1(i1)+h/2)^2 + (x2(i2)+h/2)^2 + (x3(i3)+h/2)^2); % radius of outer corner of cell
            if r1<R && r2>R % check if cell contains the source sphere surface
                for q3 = 1:length(pts) % loop over quadrature points in x3 direction
                    xq3 = x3(i3) + h*pts(q3)/2; % quadrature point x3 location (physical domain
                    for q2 = 1:length(pts) % loop over quadrature points in x2 direction
                        xq2 = x2(i2) + h*pts(q2)/2; % quadrature point x2 location (physical domain)
                        for q1 = 1:length(pts) % loop over quadrature points in x1 direction
                            xq1 = x1(i1) + h*pts(q1)/2; % quadrature point x1 location (physical domain)
                            r = sqrt(xq1^2+xq2^2+xq3^2); % quadrature point radius
                            if r <= R % check if quadrature point is inside the source sphere
                                fi = fi+1*(h/2)^3*(wts(q1)*wts(q2)*wts(q3)); % quadrature point contribution
                            end
                        end
                    end
                end
                fi = fi./(h^3); % divide integral of f over the cell by the cell volume
            elseif r2<R % check if cell completely inside source
                fi = 1;
            end
            f(i1,i2,i3) = fi; % assign cell f value
        end
    end
end
end