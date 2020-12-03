function [ y, A, At ] = sensor_operator( kData, S )
%FORWARD_MODEL Summary of this function goes here
%   Detailed explanation goes here

     % find sampling pattern (sampling mask), vectorize data
     samp = logical(abs(kData)); %[RO E1 E2 CHA PHA SET SLC]
     sampInd = find(samp~=0);
     y = kData(sampInd); % vectorized data 
     nSize = size(kData);
     A  = @(x) funA (x, sampInd, S); % image to k-space
     At = @(x) funAt(x, sampInd, S, nSize);% k-space to image

end

