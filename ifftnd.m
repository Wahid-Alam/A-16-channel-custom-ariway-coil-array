function k = ifftnd(k,fftdim,shiftflag)
% modified IFFT version specifically shift designated dimension
% default: no shift at all
% shiftdim=0, perform only fft
% specific shiftdim, shift only afterwards
% Yi Guo 04/28/2013
% modified 05/06/2013

% if fftdim==[2 3]
%     k=permute(k,[2 3 1 4 5]);
%     k=fft2(k);
%     k=permute(k,[3 1 2 4 5]);
% else
if nargin==2
for n = 1: length(fftdim)
    k=ifft(k,[],fftdim(n));
end
end

if nargin==3
for n = 1: length(fftdim)
    k=ifftshift(ifft(fftshift(k,fftdim(n)),[],fftdim(n)),fftdim(n));
end
end

% end

end
