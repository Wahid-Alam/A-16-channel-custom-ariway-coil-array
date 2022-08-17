function k = fftnd(k,fftdim,shiftflag)
% modified FFT version specifically shift designated dimension
% default: no shift at all;
% shiftdim=0, perform only fft
% specific shiftdim, shift only afterwards
% Yi Guo 04/28/2013
% modified 05/06/2013
if nargin==2
for n = 1: length(fftdim)
    k=fft(k,[],fftdim(n));
end
end

if nargin==3
for n = 1: length(fftdim)
    
    k=fftshift(fft(ifftshift(k,fftdim(n)),[],fftdim(n)),fftdim(n));
end
end
end
