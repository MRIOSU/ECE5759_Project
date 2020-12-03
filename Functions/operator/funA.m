function [y] =funA(x,sampInd,S)
% Input: x (image)
%        [RO E1 E2 CHA PHA SET SLC num_of_maps]
%        S (coil maps)
%        [RO E1 E2 CHA num_of_maps]
%        nSize (size of k-space)
% image to k-space

y = sum(bsxfun(@times, permute(S,[1 2 3 4 6 7 8 5]), x),8);%sum num_of_maps dim

for n = 1:3
  y = fftc(y,n);
end
y = y(sampInd);
