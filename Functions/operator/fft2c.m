function [x] = fft2c(x)

x = fftshift(fft(ifftshift(x,1),[],1),1)/sqrt(size(x,1));
x = fftshift(fft(ifftshift(x,2),[],2),2)/sqrt(size(x,2));