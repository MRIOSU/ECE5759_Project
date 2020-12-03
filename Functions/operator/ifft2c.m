function [x] = ifft2c(x)

x = fftshift(ifft(ifftshift(x,1),[],1),1)*sqrt(size(x,1));
x = fftshift(ifft(ifftshift(x,2),[],2),2)*sqrt(size(x,2));