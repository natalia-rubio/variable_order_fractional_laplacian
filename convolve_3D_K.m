function [I] = convolve_3D_K(K,f)
% perform convolution of K and f (0-centered)
[L1f,L2f,L3f] = size(f);
[L1K,L2K,L3K] = size(K);

K = ifftshift(K); % rearrange K so that the origin is the first element

f = padarray(f,[(L1f-1)/2 , (L2f-1)/2 , (L3f-1)/2],0,'both'); % pad f with zeros
f = ifftshift(f); % rearrange f so that the origin is the first element

FK = fftn(K)/(L1K*L2K*L3K); % 3D DFT, normalize by number of elements
FK = fftshift(FK); % rearrange FK so that the origin is at the center (for debugging)
Ff = fftn(f)./(L1f*L2f*L3f); % 3D DFT, normalize by number of elements (non-padded)
Ff = fftshift(Ff); % rearrange Ff so that the origin is at the center (for debugging)

FI = FK .* Ff; % element wise multiplication

FI_s = ifftshift(FI); % rearrange I so that the origin is the first element
I = ifftn(FI_s.*(L1K*L2K*L3K)*(L1f*L2f*L3f)); % Inverse 3D DFT, normalize by number of elements (non-padded)
I = fftshift(I); % rearrange I so that the origin is in center

I = I(L1K-3*(L1K-1)/4:L1K-(L1K-1)/4,...
    L2K-3*(L2K-1)/4:L2K-(L2K-1)/4, ...
    L3K-3*(L3K-1)/4:L3K-(L3K-1)/4); % Remove zero-padding entries
end