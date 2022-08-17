function [sigrec,ttttt] = GAspiralTTGRAPPA2(opts,kspaceUnder,acs,rseg,pseg,rep,varargin)
% SPIRALGRAPPA Through-time spiral GRAPPA reconstruction for spiral k-space
% data that has been regularly undersampled.
% 
% INPUTS
% opts  [struct]    contains info for gridding (see prepareNUFFT.m for details)
% kspaceUnder [read x proj x coils x meas] undersampled k-space
% acs [read x proj x coils x meas] fully-sampled k-space for calibrating GRAPPA weights
% rseg [1x1]    read segment size
% pseg [1x1]    projection segment size
% rep [1x1]     number of calibration measurements
% optional arguments should be entered as string-value pairs
%   'showKspace',[bool] -> set to true to display k-space during reconstruction
%   'range'             -> range for under-sampled data
% 
% OUTPUTS
% sigrec  [read x proj x coils x meas] reconstructed spiral k-space
% 
% Example: sigrec = GAspiralTTGRAPPA(opts,kspaceUnder,acs,rseg,pseg,rep)
% -----------------------------------------------------------------
% Last updated: June 2015
% MIMOSA Code Repository
% Modified by Wei-Ching Lo
% -----------------------------------------------------------------

%% Optional arguments
showKspace = false; % true-displays k-space during the GRAPPA reconstruction
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'showKspace'
            showKspace = varargin{i+1};
        case 'range'
            range = varargin{i+1};
    end
end

%% Prepare
[nr,npacc,nc,nt] = size(kspaceUnder);
np = size(acs,2);
af = np/npacc; % acceleration factor

if exist('range','var')  
    pidxacq = range;
else
    pidxacq = 1:af:np;
end

RSEG = -floor(rseg/2) : floor( (rseg-1)/2 );
PSEG = -floor(pseg/2) : floor( (pseg-1)/2 );

kernelReps=rseg*pseg*rep;
overdetFctr=kernelReps/(6*nc);
fprintf('Overdetermined factor: %.1f\n',overdetFctr);

sigrec = zeros(nr,np,nc,nt,'single');
sigrec(:,pidxacq,:,:) = kspaceUnder;
clear kspaceUnder

pidxacq(pidxacq == 0) = []; % find index for missing trajectories
pidxacq = sort(pidxacq);
pidxacq = [pidxacq pidxacq(1)+np];

sigrec = [sigrec zeros(nr,pidxacq(end)-np,nc,nt,'single')];
sigrec(:,np+1:pidxacq(end),:,:) = sigrec(:,1:pidxacq(end)-np,:,:);

acs = [acs zeros(nr,pidxacq(end)-np,nc,rep,'single')];
acs(:,np+1:pidxacq(end),:,:) = acs(:,1:pidxacq(end)-np,:,:);

% Periodically extend the trajectory, "kx" and "ky", along the projection dimension
kxbig = zeros(nr,pidxacq(end));
kybig = zeros(nr,pidxacq(end));
kxbig(:,1:np) = opts.kx;
kybig(:,1:np) = opts.ky;
kxbig(:,np+1:pidxacq(end)) = kxbig(:,1:pidxacq(end)-np);
kybig(:,np+1:pidxacq(end)) = kybig(:,1:pidxacq(end)-np);

if showKspace
    fig=figure(100); clf; set(fig,'windowstyle','docked');
    imagesc(log(makesos(sigrec(:,:,:,round(end/2)),3))); axis square
end
ttttt=[];
jj=1;
% Main reconstruction loop!
for projind = 1:length(pidxacq)-1
    proj = pidxacq(projind);

    for read = opts.gridReadLow:opts.gridReadHigh
        
        for afloop = 1:pidxacq(projind+1)-pidxacq(projind)-1%af-1
            % Location of one target point
            locX = kxbig(read,proj+afloop);
            locY = kybig(read,proj+afloop);
            % For each spiral arm, find the index of the readout point that 
            % is closest to the target point. Store this in "X".
            diffX = kxbig(opts.gridReadLow:opts.gridReadHigh,:) - locX;
            diffY = kybig(opts.gridReadLow:opts.gridReadHigh,:) - locY;
            difftot = diffX.^2 + diffY.^2;
            [~,X] = min(difftot,[],1);
            X = X + opts.gridReadLow - 1;
            
            % Calibration
            repind = 1;
            for u = 1:pseg % loops over k-space segment (projections)
                       
                psrcL = PSEG(u) + proj; %
                psrcR = pidxacq(projind+1); %
                ptarg = psrcL + afloop;
                for v = 1:rseg % loops over k-space segment (readout)
                    rsrcL = RSEG(v) + X(psrcL) - 1 : 1 : RSEG(v) + X(psrcL) + 1;
                    rsrcR = RSEG(v) + X(psrcR) - 1 : 1 : RSEG(v) + X(psrcR) + 1;
                    rtarg = RSEG(v) + X(ptarg);
                    src(:,repind:repind+rep-1) = [reshape(acs(rsrcL,psrcL,:,:),3*nc,rep); reshape(acs(rsrcR,psrcR,:,:),3*nc,rep)];
                    targ(:,repind:repind+rep-1) = reshape(acs(rtarg,ptarg,:,:),nc,rep);
                    repind = repind + rep;
                end
            end
            tic
            ws = pinv(src.')*targ.'; % GRAPPA weight set for this particular location in k-space
WS(:,:,jj)=ws;
            % Synthesis (apply GRAPPA weights to fill in missing data)
                                    ttttt=[ttttt,toc];

            psrcL = proj;
            psrcR = pidxacq(projind+1);
            ptarg = proj + afloop;
            rsrcL = X(psrcL) - 1 : 1 : X(psrcL) + 1;
            rsrcR = X(psrcR) - 1 : 1 : X(psrcR) + 1;
            rtarg = X(ptarg);
            srcr = [reshape(sigrec(rsrcL,psrcL,:,:),3*nc,nt); reshape(sigrec(rsrcR,psrcR,:,:),3*nc,nt) ];
            
            targr = (srcr.'*ws).';

            sigrec(rtarg,ptarg,:,:) = reshape(targr,[1 1 nc nt]);
            jj=jj+1;
            
        end
        if showKspace && (mod(read,300)==0)
            figure(100);imagesc(log(makesos(sigrec(:,:,:,round(end/2)),3))); axis square
            title(sprintf('Spiral GRAPPA: proj %d, read %d',proj,read),'fontname','cambria','fontsize',14,'fontweight','bold'); 
            drawnow;
        end
    end
end

sigrec = sigrec(:,[np+1:np+pidxacq(1)-1 pidxacq(1):np],:,:);

if showKspace
    imagesc(log(makesos(sigrec(:,:,:,round(end/2)),3))); axis square
    title(sprintf('Spiral GRAPPA: proj %d, read %d',proj,read),'fontname','cambria','fontsize',14,'fontweight','bold');
end

