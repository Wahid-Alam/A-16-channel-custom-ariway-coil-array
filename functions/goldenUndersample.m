function x = goldenUndersample(kspace,af)

[nr,np,nc,nt] = size(kspace);
npacc = np/af;
nttot = nt*af;
x = zeros(nr,npacc,nc,nttot,'single');
j = 0;
for i = 1:nttot
    j = j + 1;
    if j > af
        j = 1;
    end
    projIndex = (j-1)*npacc+1:j*npacc;
    x(:,:,:,i) = kspace(:,projIndex,:,ceil(i/af));
end