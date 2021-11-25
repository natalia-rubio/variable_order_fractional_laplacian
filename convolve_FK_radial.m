function [I] = convolve_radial_FK(FK_radial,f,dom)
% convert an radial FK octant to the origin-first full domain (flip over each axis)
FK_full = cat(1,FK_radial,flip(FK_radial(2:end,:,:),1));
FK_full = cat(2,FK_full,flip(FK_full(:,2:end,:),2));
FK_radial_full = cat(3,FK_full,flip(FK_full(:,:,2:end),3));
FK_radial_full = FK_radial_full./(2*dom)^3; % normalize by domain size

[L1f,L2f,L3f] = size(f); % size of f array
[L1K,L2K,L3K] = size(FK_radial_full); % size of K array

f = padarray(f,[(L1f-1)/2 , (L2f-1)/2 , (L3f-1)/2],0,'both'); % pad f with zeros
f = ifftshift(f); % rearrange f so that the origin is the first element
Ff = fftn(f)./(L1f*L2f*L3f); % 3D DFT, normalize by number of elements (non-padded)

FI_s = FK_radial_full .* Ff; % element wise multiplication
I = ifftn(FI_s.*(L1f*L2f*L3f).*((L1K*L2K*L3K))); % Inverse 3D DFT, normalize by number of elements (non-padded)
I = fftshift(I); % rearrange I so that the origin is in center
I = I(L1K-3*(L1K-1)/4:L1K-(L1K-1)/4,...
    L2K-3*(L2K-1)/4:L2K-(L2K-1)/4, ...
    L3K-3*(L3K-1)/4:L3K-(L3K-1)/4); % Remove zero-padding entries
end