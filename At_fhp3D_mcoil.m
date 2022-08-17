function At = At_fhp3D_mcoil(z, S, nx,ny,nz,nc,smaps)

z=double(z); S = double(S);

p=zeros(nx,ny,nz,nc);
p(S)=z;

for c = 1:nc
 p(:,:,:,c)=  sqrt(nx*ny*nz)*ifftnd(squeeze(p(:,:,:,c)),[1 2 3],1).*conj(squeeze(smaps(:,:,:,c)));
%p(:,:,c)= sqrt(nx*ny)*ifft2(squeeze(p(:,:,c))).*conj(squeeze(smaps(:,:,c)));
 
end

At = squeeze(sum(p,4));
