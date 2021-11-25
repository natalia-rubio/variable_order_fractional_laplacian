function [K_ramp] = K_ramp(x1 , x2 , x3)
% ramp function to suppress boundary values of K that cause instabilities
% in fft(K)
[X1,~,~] = meshgrid(x1,x2,x3);
K_ramp = 0.*X1; % Initialize K matrix
for i3 = 1:length(x3) % loop over cells in x3 direction
    for i2 = 1:length(x2) % loop over cells in x2 direction
        for i1 = 1:length(x1) % loop over cells in x1 direction
            K_ramp(i1,i2,i3) = 1-erf(4*sqrt(x1(i1)^2 + x2(i2)^2 + x3(i3)^2)/(sqrt(x1(end)^2 + x2(end)^2 + x3(end)^2)));
        end 
    end
end
end

