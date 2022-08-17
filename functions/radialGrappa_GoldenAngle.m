function [imrec,sigrec] = radialGrappa_GoldenAngle(opts,kspace,acs,rseg,pseg,rep,varargin)
% RADIALGRAPPA_GOLDENANGLE Through-time radial GRAPPA reconstruction for dynamic 2D datasets.
% Calibration data should be sorted golden angle between 0 and 180 degrees.
% Dynamic data should be unsorted golden angle between 0 and 180 degrees.
%
% INPUTS
% opts          [struct]                     contains info for gridding (see prepareNUFFT.m)
% kspace        [READ x PROJ x COILS x MEAS] undersampled k-space
% acs           [READ x PROJ x COILS x REPS] fully-sampled k-space for calibrating GRAPPA weights
% rseg          [1x1]       read segment size
% pseg          [1x1]       projection segment size
% rep          [1x1]        number of calibration frames
% optional arguments should be entered as string-value pairs
%   'projectionPadding',[string] -> 'nextFrame' or 'currentFrame'
%           controls how padding is performed along the projection dimension (default = 'currentFrame')
%   'showKspace',[bool] -> set to true to display k-space during reconstruction
%   'verbose',[bool] -> set true to show percent complete while gridding
%
% OUTPUTS
% imrec     [M x N x MEAS] reconstructed images after adaptive coil combination
% sigrec    [READ x PROJ x COIL x MEAS] reconstructed k-space
%
% EXAMPLE: imrec = radialGrappa_GoldenAngle(opts,kspace,acs,8,4,20);
% -----------------------------------------------------------------------
% Jesse Hamilton
% Dec 2013
% MIMOSA Code Repository
% -----------------------------------------------------------------------

verbose = true;
showKspace = true;
projectionPadding = 'currentFrame';
gridImage = true;

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'verbose'
            verbose = varargin{i+1};
        case 'showKspace'
            showKspace = varargin{i+1};
        case 'projectionPadding'
            projectionPadding = varargin{i+1};
        case 'gridImage'
            gridImage = varargin{i+1};
        otherwise % skip it
    end
end

[nr,npacc,nc,nttot] = size(kspace);
np = size(acs,2);
af = np/npacc;
nt = floor(nttot/af);

phi = (sqrt(5)+1)/2;

goldenAngle = 360/phi;



angles = (0:goldenAngle:(np-1)*goldenAngle)';
angles = rem(angles,360);
anglesSorted = sort(angles);




%%

kernelReps = rseg*pseg*rep;
overdetFctr=kernelReps/(6*nc);
fprintf('GRAPPA weights overdetermined by: %.1f times\n',overdetFctr);

% Pad the autocalibration matrix along projection and read dimensions
RSEG = -floor(rseg/2) : floor( (rseg-1)/2 );
PSEG = -floor(pseg/2) : floor( (pseg-1)/2 );
dr = max(abs(RSEG));

% %% Padding
% Figure out how much we padding we need along the projection dimension
angles = reshape(angles,[npacc af]);
P = zeros(size(angles));
pidxacq = zeros(size(angles));
for i = 1:af
    for j = 1:npacc
        P(j,i) = find(anglesSorted == angles(j,i));
    end
    pidxacq(:,i) = sort(P(:,i));
end

D = diff(pidxacq,1,1);
D(end+1,:) = np + pidxacq(1,:) - pidxacq(end,:);
dp = max(D(:)) + max(abs(PSEG));

angleOffset = np - pidxacq(end,:) + pidxacq(1,:);

temp = pidxacq;
pidxacq = zeros(npacc+2,af);
pidxacq(2:end-1,:) = temp;
clear temp

pidxacq(1,:) = pidxacq(2,:) - angleOffset;
pidxacq(end,:) = pidxacq(end-1,:) + angleOffset;

% %% Pad the autocalibration data
% Pad the ACS matrix along readout and projection dimensions
x = acs(:,:,:,1:rep); clear acs; % Only use the last "rep" repetitions for calibration

acs = zeros(nr+2*dr,np+2*dp,nc,rep,'double');
acs(dr+1:dr+nr,dp+1:dp+np,:,:) = x;

if strcmp(projectionPadding,'nextFrame')
    acs(dr+2:dr+nr,1:dp,:,2:end) = x(end:-1:2,end-dp+1:end,:,1:end-1);
    acs(dr+2:dr+nr,1:dp,:,1) = x(end:-1:2,end-dp+1:end,:,1);
    acs(dr+2:dr+nr,dp+np+1:end,:,1:end-1) = x(end:-1:2,1:dp,:,2:end);
    acs(dr+2:dr+nr,dp+np+1:end,:,end)=x(end:-1:2,1:dp,:,end);
elseif strcmp(projectionPadding,'currentFrame')
    acs(dr+2:dr+nr,1:dp,:,:) = x(end:-1:2,end-dp+1:end,:,:);
    acs(dr+2:dr+nr,dp+np+1:end,:,:) = x(end:-1:2,1:dp,:,:);
else error('Projection padding method not defined correctly')
end

clear x
acs(1:dr,:,:,:) = acs(2*dr:-1:1+dr,:,:,:);
acs(dr+nr+1:end,:,:,:) = acs(nr+dr:-1:nr+1,:,:,:);

%% Initialize matrix to hold reconstructed k-space
sigrec = zeros(nr,np+2*dp,nc,nttot,'double');
for u = 1:af
    for p = 1:npacc
        sigrec(:,dp+P(p,u),:,u:af:end) = kspace(:,p,:,u:af:end);
    end
    sigrec(2:end,1:dp,:,u:af:end) = sigrec(end:-1:2,np+1:np+dp,:,u:af:end);
    sigrec(2:end,dp+np+1:end,:,u:af:end) = sigrec(end:-1:2,dp+1:2*dp,:,u:af:end);
end

%% Reconstruct
anglegaps = diff(pidxacq);
gapsizes = unique(anglegaps(:));
prevlen = 0;
if showKspace; fig=figure(100); clf; end
for ig = 1:length(gapsizes)
    accel = gapsizes(ig);
    
    read = repmat([1 3 5 1 3 5],1,nc)';
    read = repmat(read,1,rseg);
    read = bsxfun(@plus,read,RSEG);
    read = read(:);
    read = repmat(read,1,pseg);
    read = read(:);
    read = repmat(read,rep,1);
    read = read - RSEG(1);
    
    proj = repmat([1 1 1 accel+1 accel+1 accel+1],1,nc)';
    proj = repmat(proj,1,rseg);
    proj = reshape(proj,[],1);
    proj = repmat(proj,1,pseg);
    proj = bsxfun(@plus,proj,PSEG);
    proj = proj(:);
    proj = repmat(proj,rep,1);
    proj = proj - PSEG(1);
    
    coil = reshape( repmat(1:nc,6,1), 6*nc, 1);
    coil = repmat(coil,1,pseg*rseg*rep);
    coil = coil(:);
    
    meas = repmat(1:rep,6*nc*rseg*pseg,1);
    meas = meas(:);
    
    matsize = [rseg+4,accel+pseg,nc,rep];
    
    acsIdx = sub2ind(matsize,read,proj,coil,meas);
    
    read = repmat([1 3 5 1 3 5],1,nc)';
    read = repmat(read,nt,1);
    proj = repmat([1 1 1 accel+1 accel+1 accel+1],1,nc)';
    proj = repmat(proj,nt,1);
    coil = reshape( repmat(1:nc,6,1), 6*nc, 1);
    coil = repmat(coil,nt,1);
    meas = repmat(1:nt,6*nc,1);
    meas = meas(:);
    
    matsize = [5,accel+1,nc,nt];
    
    sigIdx = sub2ind(matsize,read,proj,coil,meas);
    
    % Target points
    read = repmat(ones(1,accel-1),1,nc)';
    read = repmat(read,1,rseg);
    read = bsxfun(@plus,read,RSEG);
    read = read(:);
    read = repmat(read,1,pseg);
    read = read(:);
    read = repmat(read,rep,1);
    read = read-RSEG(1);
    
    proj = repmat(2:accel,1,nc)';
    proj = repmat(proj,1,rseg);
    proj = reshape(proj,[],1);
    proj = repmat(proj,1,pseg);
    proj = bsxfun(@plus,proj,PSEG);
    proj = proj(:);
    proj = repmat(proj,rep,1);
    proj = proj - PSEG(1);
    
    coil = reshape( repmat(1:nc,(accel-1),1), (accel-1)*nc, 1);
    coil = repmat(coil,1,pseg*rseg*rep);
    coil = coil(:);
    
    meas = repmat(1:rep,(accel-1)*nc*rseg*pseg,1);
    meas = meas(:);
    
    matsize = [rseg+4,accel+pseg,nc,rep];
    
    acsIdxTarg = sub2ind(matsize,read,proj,coil,meas);
    
    read = repmat(ones(1,accel-1),1,nc)';
    read = repmat(read,nt,1);
    proj = repmat(2:accel,1,nc)';
    proj = repmat(proj,nt,1);
    coil = reshape( repmat(1:nc,accel-1,1), (accel-1)*nc, 1);
    coil = repmat(coil,nt,1);
    meas = repmat(1:nt,(accel-1)*nc,1);
    meas = meas(:);
    
    matsize = [5,accel+1,nc,nt];
    
    sigIdxTarg = sub2ind(matsize,read,proj,coil,meas);
    clear coil proj read meas matsize
    
    for u=1:af
        pidx = find(anglegaps(:,u)==accel);
        for itp = 1:length(pidx)
            p = pidxacq(pidx(itp),u);
            msg = sprintf('Golden Angle Radial GRAPPA: interleaf %d/%d, projection %d/%d\n',u,af,p,np);
            fprintf([repmat('\b',1,prevlen) '%s'],msg);
            prevlen = numel(msg);
            for r = 3:nr-2
                calMatrix = acs(r-2+RSEG(1)+dr:r+2+RSEG(end)+dr,p+PSEG(1)+dp:p+dp+PSEG(end)+accel,:,:);
                src = reshape(calMatrix(acsIdx),6*nc,rseg*pseg*rep);
                targ = reshape(calMatrix(acsIdxTarg),(accel-1)*nc,rseg*pseg*rep);
                ws = pinv(src.')*targ.';
                
                recMatrix = sigrec(r-2:r+2,dp+p:dp+p+accel,:,u:af:end);
                srcr = reshape(recMatrix(sigIdx),6*nc,nt);
                recMatrix(sigIdxTarg) = reshape((srcr.'*ws).',[1 accel-1 nc nt]);
                sigrec(r-2:r+2,dp+p:dp+p+accel,:,u:af:end) = recMatrix;
            end
            if showKspace
                imagesc(log(makesos(sigrec(:,:,:,u+af),3))); axis square; axis off;
                title(sprintf('Interleaf %d',u),'fontsize',12,'fontweight','bold','fontname','cambria'); drawnow;
            end
        end
    end
end
sigrec = sigrec(:,dp+1:dp+np,:,:);
clear acs

if showKspace, close(fig); end

%% Grid
imrec = [];
if gridImage
    if ~strcmp(opts.viewOrder,'goldenAngle_sorted_180')
        opts = prepareNUFFT(opts.N,np,'radial','goldenAngle_sorted_180');
    end
    imrec = kspaceToImage(opts,sigrec,'verbose',verbose);
end