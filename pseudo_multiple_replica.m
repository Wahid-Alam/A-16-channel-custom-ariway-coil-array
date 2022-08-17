function [g_pseudo, snr_pseudo]=pseudo_multiple_replica(data, smask,smap)
% smask = zeros(nx,ny); 
% smask(1:R:end,:)=1;
% smask = fftshift(smask);
% smask = repmat(smask,[1 1 nc]);
[nx ny nz nc]= size(smask);
S=double(find(logical(smask)));
A = @(z)A_fhp3D_mcoil(z, S,nx,ny,nz,nc,smap);
At = @(z)At_fhp3D_mcoil(z, S, nx,ny,nz,nc,smap);
opts.nx = nx; 
opts.ny = ny; 
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
img_noise_rep = zeros(nx,ny,nz,50); 
for r=1:50
    r
    noise_white = complex(randn(size(smask)),randn(size(smask)));
    noise_white = noise_white.*smask;
    s = data + double(noise_white(S));
    x_init = At(s); 
    grad_right = At(s);grad_right = grad_right(:);
    x0=x_init(:);
    %x0=rand(size(x0));
    [x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,20,[],[],x0(:));
    x_recon=reshape(x_recon,opts.nx, opts.ny,opts.nz);
    img_noise_rep(:,:,:,r) = x_recon;
end

g_pseudo = std(abs(img_noise_rep + max(abs(img_noise_rep(:)))),[],4);
g_pseudo(g_pseudo < eps) = 1;
snr_pseudo = mean(img_noise_rep,4)./g_pseudo;
g_pseudo = g_pseudo.*sqrt(sum(abs(smap).^2,4));
