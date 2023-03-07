function Y = remove_RO_oversamling( data )
%   data: firt dimension should be Readout dimension

fd(1) = 0.25; fd(2) = 0.25;
N=size(data);
I = fftshift(ifft(ifftshift(data,1),[],1),1)* sqrt(N(1));
I(1:round(N(1)*fd(1)),:,:,:,:,:,:,:,:,:)=[]; % remove first 1/4
I(end-round(N(1)*fd(2))+1:end,:,:,:,:,:,:,:,:,:)=[]; % remove last 1/4
Y = fftshift(fft(ifftshift(I,1),[],1),1)/sqrt(size(I,1));% / sqrt(N(1)/(N(1)-z1-z2));
RO = size(data,1);
samp = logical(sum(abs(reshape(data,[RO, numel(data)/RO])),2));
pre_z = find(samp == 0, 1, 'last' ); % parameter related to assymetry echo
pre_z = floor(pre_z*(1-sum(fd)));
Y(1:pre_z,:,:,:,:,:,:,:,:) = 0; %make sure the unsampled data zero

end

