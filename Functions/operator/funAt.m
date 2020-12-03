function [y] =funAt(x,sampInd, S, nSize)
% Input: x (vectorized k-space)
%        S (coil maps)
%        [RO E1 E2 CHA num_of_maps]
%        nSize (size of k-space, RO E1 E2 CHA PHA SET SLC)
% k-space to image
y = zeros(nSize);
y(sampInd) = x;
for n = 1:3
  y = ifftc(y,n);
end
y =  sum(bsxfun(@times, permute(conj(S),[1 2 3 4 6 7 8 5]), y), 4); %sum coil dim
