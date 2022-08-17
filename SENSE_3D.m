%% extraction of k-space from P-file
%clear all
%breathing - gre
raw_data_file = '/Users/moalam/Desktop/Research_folder/basic_parallelMRI/3D_pseudo/Dec10_scans/P45568.7';
addpath('/Users/moalam/Desktop/Research_folder/basic_parallelMRI/orchestra-sdk-1.10-1.matlab/');

%addpath(fullfile(getenv('TOOLBOX_PATH'), 'matlab')); % add BART tool bax path


%read pertinent pfile header data
pfile = GERecon('Pfile.Load', raw_data_file);
header = GERecon('Pfile.Header', pfile);
acquiredSlices = pfile.slicesPerPass;
outputSlices = pfile.reconstructedSlicesPerPass;
scaleFactor = outputSlices / pfile.scaleFactor3d;
SIZE_Z_MATRIX = acquiredSlices;
SIZE_Y_MATRIX = header.RawHeader.user2;
SIZE_X_MATRIX = pfile.xRes;

%read acquisition order from external file.
%acq_order = int32(dlmread(capr_view_order)); %sjk better way to do this?

%error condition checks
if (pfile.echoes ~= 1)
    error('mum echoes must be 1')
end
if (pfile.passes ~= 1)
    error('num passes must be 1')
end
echo = 1;
pass = 1;

%load all raw data
for slice = 1:acquiredSlices
    
    sliceInfo.pass = pass;
    sliceInfo.sliceInPass = slice;
    
    for channel = 1:pfile.channels
        
        % Load K-Space
        kSpace(:,:,slice,channel) = GERecon('Pfile.KSpace', sliceInfo, echo, channel);
        
    end
end

[nx,ny,nz,nc]=size(kSpace);

%%
kSpace= k_out;  
[nx,ny,nz,nc]=size(kSpace);
mc_imgs = ifftnd(kSpace,[1 2 3],1);                                                                  
                                    
                                    
% module to create low resoultion images for coil map estimation

sigmaf=4;%sigma for a Gaussian mask to perform low pass filtering n k-space; adjust accordingly. 

%Roughly sigmaf=4 corresponds to using 5% of the center of k-space
% lowres_mask= fspecial('gaussian', [nx ny], sigmaf);
% lowres_mask=lowres_mask./max(lowres_mask(:));
 
kspace_decoupled =  ifftnd(kSpace,[3],1); %bart('fft -u -i 1', kSpace);

lowres_mask=create_lpf(nx,ny,sigmaf);


% for i = 1:pfile.channels,
for i = 1:nc
    for j = 1:nz
        kspace_lowres_decoupled(:,:,j,i)= squeeze(kspace_decoupled(:,:,j,i)).*lowres_mask;
    end
end

% create low resolution images for csm estimation

mc_lowresimgs = ifftnd(kspace_lowres_decoupled,[1 2],1); 

% Sum of squares CSM map estimation
im_rss_lowres = sqrt(sum(abs(mc_lowresimgs).^2,4));

for i = 1:nc
sens_maps(:,:,:,i) = (mc_lowresimgs(:,:,:,i))./(im_rss_lowres);
end
%sens_maps=ifftshift(sens_maps,3);
%
%%
Rx=1;Ry=1; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;

S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
 
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);

% First guess, direct IFFT
x_init = At(b);

% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:);
 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_R_1x1x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
%x_recon_h_side_R_1x1x1= ifftshift(x_recon_h_side_R_1x1x1, 3);
% % generating SNR and g maps from the pseudo replica method
% [g_R_1x1x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

%%
figure();
subplot(1,3,1); imagesc(abs((squeeze(snr_h_R5(:,:,end)))));colorbar; title('after prewhitening SNR map at R=5');axis image;
subplot(1,3,2); imagesc(abs(squeeze(g_h_R5(:,:,end)))); colormap(jet); colorbar; title('after prewhitening g-map at R=5');axis image;
subplot(1,3,3); imagesc(abs((squeeze( x_recon_h_R5(:,:,end)))));colorbar; title('3D SENSE recon R=5');axis image;
%%
figure();
subplot(1,2,1); imagesc(abs(squeeze(g_h_prew_R4(:,:,end)))); colormap(jet); colorbar; title('after prewhitening g-map at R=5');axis image;
subplot(1,2,2); imagesc(abs((squeeze(x_recon_h_R4(:,:,end)))),[0 50]);colorbar; title('3D SENSE recon R=5');axis image;
%%
figure();
imagesc(abs(squeeze(x_recon_aw_R_1x1x1(:,:,end))));colormap(gray);
head_mask= roipoly();
%%
figure();
imagesc(abs(squeeze(x_recon_aw_R_1x1x1(:,:,end))));colormap(gray);
bcgrnd_mask= roipoly();
%%
for i=1:16
    masked= squeeze(mc_imgs(:,:,end,i)).* bcgrnd_mask;
    prew_wts(i)=var(masked(:));
end

%%
indv_coil_maps= [sens_maps(:,:,end,1).*head_mask sens_maps(:,:,end,2).*head_mask sens_maps(:,:,end,3).*head_mask sens_maps(:,:,end,4).*head_mask sens_maps(:,:,end,5).*head_mask sens_maps(:,:,end,6).*head_mask sens_maps(:,:,end,7).*head_mask sens_maps(:,:,end,8).*head_mask...
    sens_maps(:,:,end,9).*head_mask sens_maps(:,:,end,10).*head_mask sens_maps(:,:,end,11).*head_mask sens_maps(:,:,end,12).*head_mask sens_maps(:,:,end,13).*head_mask sens_maps(:,:,end,14).*head_mask sens_maps(:,:,end,15).*head_mask sens_maps(:,:,end,16).*head_mask ];

ismrm_imshow(angle(indv_coil_maps)); brighten(0.4); colormap(gray); colorbar;
%%
indv_coil_imgs= [mc_imgs(:,:,end,1)/prew_wts(1) mc_imgs(:,:,end,2)/prew_wts(2) mc_imgs(:,:,end,3)/prew_wts(3) mc_imgs(:,:,end,4)/prew_wts(4) mc_imgs(:,:,end,5)/prew_wts(5) mc_imgs(:,:,end,6)/prew_wts(6) mc_imgs(:,:,end,7)/prew_wts(7) mc_imgs(:,:,end,8)/prew_wts(8)...
    mc_imgs(:,:,end,9)/prew_wts(9) mc_imgs(:,:,end,10)/prew_wts(10) mc_imgs(:,:,end,11)/prew_wts(11) mc_imgs(:,:,end,12)/prew_wts(12) mc_imgs(:,:,end,13)/prew_wts(13) mc_imgs(:,:,end,14)/prew_wts(14) mc_imgs(:,:,end,15)/prew_wts(15) mc_imgs(:,:,end,16)/prew_wts(16) ];

ismrm_imshow(abs(indv_coil_imgs));brighten(0.2); colormap(gray); colorbar;
%%

S =double(find(logical(mask_4)));
b=double(k_out(S));

A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_R_2x2x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
% [g_aw_R_2x2x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);
clear b;

S =double(find(logical(mask_6)));
b=double(k_out(S));

A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_R_3x2x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
% [g_aw_R_2x2x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);
clear b;

S =double(find(logical(mask_8)));
b=double(k_out(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_R_4x2x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
% [g_aw_R_2x2x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);
clear b;
%%
Rx=3;Ry=1; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_3x1x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_3x1x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

Rx=4;Ry=1; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_4x1x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_4x1x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

Rx=5;Ry=1; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_5x1x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_5x1x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

Rx=1;Ry=2; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_1x2x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_1x2x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);


Rx=1;Ry=3; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_1x3x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_1x3x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

Rx=1;Ry=4; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_1x4x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_1x4x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

Rx=1;Ry=5; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_1x5x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_1x5x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);


Rx=2;Ry=2; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_2x2x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_2x2x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

Rx=2;Ry=3; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_2x3x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_2x3x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

Rx=3;Ry=2; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_3x2x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_3x2x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

Rx=4;Ry=2; Rz=1; 
sampling_mask = zeros(size(kSpace));  
sampling_mask(1:Rx:end,1:Ry:end,1:Rz:end,:)=1;
S =double(find(logical(sampling_mask)));
b=double(kSpace(S));
A = @(z)A_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
At = @(z)At_fhp3D_mcoil(z,S,nx,ny,nz,nc,sens_maps);
opts.nx = nx; 
opts.ny = ny;
opts.nz = nz;
opts.nc = nc; 
grad_left = @(z)Eval_grad_left(z,A,At,opts);
% First guess, direct IFFT
x_init = At(b);
% write code for simple SENSE recon
grad_right = At(b);grad_right = grad_right(:); 
x0=x_init(:);
[x_recon,~,~,~,resvec] = pcg(grad_left,grad_right(:),1e-5,30,[],[],x0(:));
x_recon_aw_R_4x2x1=reshape(x_recon,opts.nx, opts.ny,opts.nz);
% % generating SNR and g maps from the pseudo replica method
[g_aw_R_4x2x1, ~]=pseudo_multiple_replica(b, sampling_mask,sens_maps);

