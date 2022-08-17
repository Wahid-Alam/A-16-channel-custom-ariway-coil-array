%% example script to read pfiles from the premier. 
% note - this is a quick and dirty script; - Sajan 

%%
addpath('nufft_toolbox/');
addpath('fidall_orig/');

%% specify the pfile and the fidall waveform location

NARMS=34;TR =5.0e-3; 
if nav==0
reconfile = ('/Users/moalam/Desktop/Research_folder/grappa_through_time/grappa_code/waveforms/210_83_34_UDS_gr_IMAG_5400_smax_130_rofi'); %wahid: you need to change this to vds_21_rofi
elseif nav==1
    reconfile = ('/Users/moalam/Desktop/Research_folder/grappa_through_time/grappa_code/waveforms/210_83_34_UDS_gr_NAV_5400_smax_130_rofi'); %wahid: you need to change this to vds_21_rofi
  end
%% read the Pfile using Rolfs scripts.

[dd,k,wf]=recon_spiral_stripped(pfile, 'nix',reconfile,0,0,0,false,pfile,false);

dd = dd(:,:,:);

%% some data reshaping

dd=dd(:,:);
mtx_reco =[84,84]; 
[nviews,nr]=size(wf.t); 
[nc, ~]=size(dd); 

dd = reshape(dd, nc, nr, nviews); 
k = reshape(k, nr, nviews); 
dcf = reshape(wf.dcf, nr, nviews); 


nc=size(dd,1); 

% data dimensions
kdata = permute(dd,[2 3 1]);
kdata=kdata(:,1:end,:);
[nx,ntviews,nc]=size(kdata);
kloc = k;
w=dcf;
load dcf_630readouts.mat;

% Image matrix size
n1=84; 
n2= 84; 

% no of desired spiral arms       
narms =34;

% number of frames
nt=floor(ntviews/narms);


% crop the data according to the number of arms per frame
kdata=kdata(:,1:nt*narms,:);
kdata = 1.4790e3*kdata./max(abs(kdata(:)));
k=kloc(:,1:nt*narms);
w=w(1:size(k,1),1);
w=repmat(w,[1,nt*narms]);