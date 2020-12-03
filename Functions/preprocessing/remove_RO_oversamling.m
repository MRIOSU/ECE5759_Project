function Y = remove_RO_oversamling( data )
%   data: firt dimension should be Readout dimension

fd(1) = 0.25; fd(2) = 0.25;
N=size(data);
I = fftshift(ifft(ifftshift(data,1),[],1),1)* sqrt(N(1));
I(1:round(N(1)*fd(1)),:,:,:,:,:,:,:,:,:)=[]; % remove first 1/4
I(end-round(N(1)*fd(2))+1:end,:,:,:,:,:,:,:,:,:)=[]; % remove last 1/4
Y = fftshift(fft(ifftshift(I,1),[],1),1)/sqrt(size(I,1));% / sqrt(N(1)/(N(1)-z1-z2));

end

