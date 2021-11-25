T = 0.1;
Fs = 1/T;
D = 10;
x = -D:T:D;
x1 = 0:T:D;
[X1,X2,X3] = meshgrid(x,x,x);
L = length(x);
cen = ceil(L/2);

r = sqrt(3*((0:T:D).^2));
L_r = length(r);

%s = @(r) cos(x*2*pi);
%s = @(r) 0.*r + 1;
%s = @(r) 1./(sqrt(2*pi)) .* exp(-0.5*r.^2);
figure(3)
s3d = @(x1,x2,x3) 100./(sqrt(2*pi)) .* exp(-0.7*sqrt(x1.^2 + x2.^2 + x3.^2).^2);
K = ifftshift(s3d(X1,X2,X3));
%surf(X1(:,:,cen),X2(:,:,cen),K(:,:,1))
F = fftn(ifftshift(s3d(X1,X2,X3)));
f = (Fs/(L)).*((0:L/2));
f = [f,-flip(f(2:end))];
figure(1)
%plot(f, F(:,1,1).* ((2*D)/L)^3*2*pi/3) ; %.*  ((4*pi*D^3/3)/((2*D)^3)))
%plot(f, F(:,1,1).* ((4*pi*(max(r))^3/3)/(L^3))) %.* ((4*pi*(max(r))^3/3)/((2*D)^3)))
F_norm = F .* ((2*D)^3)/(L^3);
plot(f, F_norm(:,1,1)) %.* ((4*pi*(max(r))^3/3)/((2*D)^3)))
%%
disp("FK_3D_norm")
norm(F_norm(:))
%figure(4)
%%


%s = @(r) cos(r*2*pi);
%s = @(r) 0.*r + 1;
s = @(r) 100./(sqrt(2*pi)) .* exp(-0.7*r.^2);
T_int_rad = max(r);
[FK_radial, k] = compute_FK_radial_uni(s,r,T,T_int_rad) ;
FK_radial = FK_radial;
%k = [k, -flip(k(2:end))];
%FK_radial = [FK_radial];
figure(2)
plot(k,FK_radial)
[FK_radial_3D, K_radial_3D] = map_radial_to_3D(FK_radial, T, r, x1, x1, x1,s); % Map radial FK into 3D

FK_full = cat(1,FK_radial_3D,flip(FK_radial_3D(2:end,:,:),1));
FK_full = cat(2,FK_full,flip(FK_full(:,2:end,:),2));
FK_radial_full = cat(3,FK_full,flip(FK_full(:,:,2:end),3));

disp("FK_radial_norm")
    norm((FK_radial_full(:)))
    
    %%
    figure(3)
 surf(X1(:,:,cen),X2(:,:,cen),FK_radial_full(:,:,1))
 %surf(X1(:,:,cen),X2(:,:,cen),abs(F_norm(:,:,1)))
 figure(4)
 plot(k,FK_radial); hold on;
 plot(f, F_norm(:,1,1))
 figure(5)
 surf(X1(:,:,cen),X2(:,:,cen),FK_radial_full(:,:,1)-abs(F_norm(:,:,1)))