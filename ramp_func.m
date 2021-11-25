function [ramp] = ramp_func(u)
% Ramp function for computing k1 and k2
ramp = -(u-1)^2 + 1;
if u > 1
    ramp = 1;
end
end 