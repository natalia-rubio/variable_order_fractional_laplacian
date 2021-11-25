function [I_exp] = explicit_sol(x1,x2,x3, R, s)
% explicit solution for gravitational potential
G = 1; % gravitational constant
M = 1; % total mass
l = 1; % length constant
H = 3*gamma(3/2-s)*G*M*l^(2-2*s) / (4^(s-1/2)*sqrt(pi)*(2*s-1)*gamma(s)*(R^3)*(2*s)*(2*s+1)); % constant coefficient

[I_exp,~,~] = meshgrid(x1,x2,x3);
I_exp = 0.*I_exp;

for i3 = 1:length(x3) % loop over elements in x3 direction
    for i2 = 1:length(x2) % loop over elements in x2 direction
        for i1 = 1:length(x1) % loop over elements in x1 direction
            r = sqrt(x1(i1)^2+x2(i2)^2+x3(i3)^2); % find radial distance r
            % piecewise explicit solution
            if r == 0
                I_exp(i1,i2,i3) = 0;
            elseif  r > R
                I_exp(i1,i2,i3) = -H*r^(2*s)*...
                    ((1-R/r)^(2*s) * (1+2*s*R/r) - (1+R/r)^(2*s) * (1-2*s*R/r)); % exterior solution
            elseif r <= R
                I_exp(i1,i2,i3) = H*r^(2*s)*...
                    ((1+R/r)^(2*s) * (1-2*s*R/r) + (R/r-1)^(2*s) * (1+2*s*R/r)); % interior solution
            end
        end
    end
end
end