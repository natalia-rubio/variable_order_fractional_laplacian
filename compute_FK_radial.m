function [FK_radial, k] = compute_FK_radial(s,rK,h,domain_size)
max_iter = 20000; % set maximum of iterations
dk = 1/(length(rK) * h); % discretization frequency
k = dk*(0:(length(rK)/2)); % Fourier space discretization
FK_radial = k.*0; % initialize radial Fourier space K
lambda = 0.001; % regularization for radial Fourier transform
for i = 1:length(k)
    % compute radial Fourier space K for each frequency
    [FK_radial(i), ~] = rad_fourier_int_exp(k(i), lambda, s, max_iter, domain_size);
end
%FK_radial = [FK_radial(1), FK_radial(2:end), flip(FK_radial(2:end))];
end

