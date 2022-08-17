function imadapt = kspaceToImage(opts,kspace,varargin)
% KSPACETOIMAGE Grids multi- or single-channel non-Cartersian k-space to
% the image domain. Then runs adaptive coil combine.
%
% INPUTS
% opts [struct] contains information needed for gridding (see prepareNUFFT.m)
% kspace [READ X PROJ],[READ x PROJ x COILS],[READ x PROJ x COILS x MEAS],[READ x PROJ x COILS x PART x MEAS]
%       single or multi-channel images
% optional arguments should be input as string-value pairs
%       'verbose',[bool]  -> set to true to display percent complete (default = true)
%
% OUTPUT
% imadapt [M x N], [M x N x MEAS], [M x N x PART x MEAS]
%   Images after adaptive coil combination
%
% EXAMPLES: imadapt = kspaceToImage(opts,kspace);
% -----------------------------------------------------------------------------------------
% Jesse Hamilton
% Dec 2013
% MIMOSA Code Repository
% -----------------------------------------------------------------------------------------

imadapt = [];
verbose = 1;
bShowImages = false;
strMethod = 'Adapt';
invFFT = false;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'verbose'
            verbose = varargin{i+1};
        case 'showImages'
            bShowImages = varargin{i+1};
        case 'method'
            strMethod = varargin{i+1};
            if ~(strcmp(strMethod,'SOS') || strcmp(strMethod,'Adapt'))
                error('Invalid coil combination method');
            end
        case 'invFFT'
            invFFT = varargin{i+1};
        otherwise % skip it
    end
end

rlow = opts.gridReadLow;
rhigh = opts.gridReadHigh;
nr = rhigh-rlow+1;

np = size(opts.kx,2); % projections in fully-sampled data
npacc = size(kspace,2);
af = np/npacc
mask = true(opts.N,opts.N);

switch ndims(kspace)
    case 2
        kspace = kspace(opts.gridReadLow:opts.gridReadHigh,:);
        x = kspace(:);
        imadapt = embed(opts.G'*(opts.wib.*x(:)),mask);
        imadapt = af*imadapt;
        
    case 3
        kspace = kspace(opts.gridReadLow:opts.gridReadHigh,:,:);
        nc = size(kspace,3);
        imcoil = zeros(opts.N,opts.N,nc,'single');
        x = zeros(nr,np);
        for c = 1:nc
            x(:,1:af:end) = kspace(:,:,c);
            imcoil(:,:,c) = embed(opts.G'*(opts.wib.*x(:)),mask);
        end
        if strcmp(strMethod,'Adapt')
            imadapt = af*openadapt_wrapper(imcoil);
        elseif strcmp(strMethod,'SOS')
            imadapt = af*makesos(imcoil,3);
        else
            error('Invalid coil combination method');
        end
        
    case 4
        kspace = kspace(opts.gridReadLow:opts.gridReadHigh,:,:,:);
        [~,~,nc,nt] = size(kspace);
        imcoil = zeros(opts.N,opts.N,nc,nt,'single');
        prevlen = 0;
        x = zeros(nr,np);
        
        if sum(strcmp(opts.viewOrder,{'linear_sorted','goldenAngle_sorted_180','goldenAngle_sorted_360'}))
            for c = 1:nc
                if verbose
                    msg = sprintf('NUFFT (k-space to image): Coil %d/%d\n',c,nc);
                    fprintf([repmat('\b',1,prevlen) '%s'],msg); prevlen = numel(msg);
                end
                for t = 1:nt
                    x(:,1:af:end) = kspace(:,:,c,t);
                    imcoil(:,:,c,t) = embed(opts.G'*(opts.wib.*x(:)),mask);
                end
            end
            
        elseif sum(strcmp(opts.viewOrder,{'goldenAngle_180','goldenAngle_360','interleaved'}))
            for u = 1:af
                projIndex = (u-1)*npacc+1:u*npacc;
                x = zeros(nr,np);
                for c = 1:nc
                    if verbose
                        msg = sprintf('NUFFT (k-space to image): Interleaf %d/%d, Coil %d/%d\n',u,af,c,nc);
                        fprintf([repmat('\b',1,prevlen) '%s'],msg); prevlen = numel(msg);
                    end
                    for t = u:af:nt
                        x(:,projIndex) = kspace(:,:,c,t);
                        imcoil(:,:,c,t) = embed(opts.G'*(opts.wib.*x(:)),mask);
                    end
                end
            end
        end
        imadapt = af*openadapt_wrapper(imcoil);
        
    case 5
        kspace = squeeze(kspace(opts.gridReadLow:opts.gridReadHigh,:,:,:,:));
        [~,~,nc,nz,nt] = size(kspace);
        if invFFT
            for t = 1:nt
                fprintf('Inverse FFT along z: meas %d\n',t);
                kspace(:,:,:,:,t) = ifftshift(ifft(ifftshift(kspace(:,:,:,:,t),4),[],4),4);
            end
        end
        imcoil = zeros(opts.N,opts.N,nc,nt,'single');
        imadapt = zeros(opts.N,opts.N,nz,nt,'single');
        prevlen = 0;
        x = zeros(nr,np);
        if bShowImages, fig=figure; end
        if sum(strcmp(opts.viewOrder,{'linear_sorted','goldenAngle_sorted_180','goldenAngle_sorted_360'}))
            for itPart = 1:nz
                x = zeros(nr,np);
                for itCoil = 1:nc
                    msg = sprintf('NUFFT: Part %d/%d, Coil %d/%d\n',itPart,nz,itCoil,nc);
                    fprintf([repmat('\b',1,prevlen) '%s'],msg); prevlen = numel(msg);
                    for itMeas = 1:nt
                        x(:,1:af:end) = squeeze(kspace(:,:,itCoil,itPart,itMeas));
                        imcoil(:,:,itCoil,itMeas) = embed(opts.G'*(opts.wib.*x(:)),mask);
                    end
                end
                if strcmp(strMethod,'Adapt')
                    imadapt(:,:,itPart,:) = af*openadapt_wrapper(imcoil);
                elseif strcmp(strMethod,'SOS')
                    imadapt(:,:,itPart,:) = af*makesos(imcoil,3);
                else
                    error('Method not defined');
                end
                if bShowImages
                    imshow(squeeze(abs(imadapt(:,:,itPart,round(end/2)))),[]);
                    title(sprintf('Partition %d',itPart)); drawnow;
                end
            end
        elseif sum(strcmp(opts.viewOrder,{'goldenAngle_180','goldenAngle_360','interleaved'}))
            for z = 1:nz
                for u = 1:af
                    if verbose
                        msg = sprintf('NUFFT (k-space to image): Par %d/%d, Interleaf %d/%d\n',z,nz,u,af);
                        fprintf([repmat('\b',1,prevlen) '%s'],msg); prevlen = numel(msg);
                    end
                    projIndex = (u-1)*npacc+1:u*npacc;
                    x = zeros(nr,np);
                    for c = 1:nc
                        for t = u:af:nt
                            x(:,projIndex) = kspace(:,:,c,z,t);
                            imcoil(:,:,c,t) = embed(opts.G'*(opts.wib.*x(:)),mask);
                        end
                    end
                end
                if strcmp(strMethod,'Adapt')
                    imadapt(:,:,z,:) = af*openadapt_wrapper(imcoil);
                elseif strcmp(strMethod,'SOS')
                    imadapt(:,:,z,:) = af*makesos(imcoil,3);
                else
                    error('Method not defined');
                end
                if bShowImages
                    imshow(squeeze(abs(imadapt(:,:,z,ceil(end/2)))),[]);
                    title(sprintf('Partition %d',z)); drawnow;
                end
            end
        end
        if bShowImages
            close(fig);
        end
        
    otherwise
        error('K-space has wrong dimensions')
end

