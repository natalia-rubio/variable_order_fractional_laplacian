function [K_full] = oct2full(K_oct)
% convert an octant to the full domain (flip over each axis)
K_full = cat(1,flip(K_oct(2:end,:,:),1),K_oct);
K_full = cat(2,flip(K_full(:,2:end,:,:),2),K_full);
K_full = cat(3,flip(K_full(:,:,2:end),3),K_full);
end