function A = A_fhp3D_mcoil(z, S,nx,ny,nz,nc,smaps)
z=double(z); S = double(S);
p = zeros(nx,ny,nz,nc); 
for i = 1: nc 
    p(:,:,:,i)= 1/sqrt(nx*ny*nz) * fftnd(z.*smaps(:,:,:,i),[1 2 3],1);
    %p(:,:,i)=1/sqrt(nx*ny)*fft2(z.*smaps(:,:,i));
end
A = p(S);
end

