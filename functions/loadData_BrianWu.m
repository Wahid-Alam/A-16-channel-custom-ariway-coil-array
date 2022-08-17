function [kspace,kx,ky] = loadData_BrianWu(loadname,N,np)

% load(fullfile(datadir,'sleep48128'));
load(loadname)
kspace = k_rad; clear k_rad; % I like to call my dynamic data "kspace" :)

[nr,nptot,nc] = size(kspace);
kx = real(T(:,1:np))*N;
ky = imag(T(:,1:np))*N;
clear T
% figure; plot(kx,ky,'g'); axis image; title('Radial Traj');

kspace = permute(kspace,[1 3 2]);

nt = floor(nptot/np);
nptot = nt*np; % ensure that total number of projections is divisible by 144
kspace = kspace(:,:,end-nptot+1:end);
kspace = reshape(kspace,[nr nc np nt]);
kspace = permute(kspace,[1 3 2 4]);
fprintf('Read = %d, Proj = %d, Coils = %d, Meas = %d\n',nr,np,nc,nt);