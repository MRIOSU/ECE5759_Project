function Y = coil_compression(data, ncoil)
%input: data [RO E1 E2 CHA PHA SET SLC]
%       ncoil number of virtual coils

N=size(data);
% N(chDim),
Nc = ncoil;
chDim = 4; % chanel dimension
slDim = 7; % Slice dimension

dataC = dimData(data,[],chDim, [1:min(Nc,N(chDim))],'r'); % Initialization

nPer = 1:max(ndims(data),slDim);nPer(chDim)=nPer(end); nPer(end)=chDim;
if Nc < N(chDim) % only compress if virtual channels are lesser than actual channels
    for i=1:size(data,slDim)
        tmp = dimData(data,[],slDim,i,'r'); % Read data along 'slDim'
        tmp = permute(tmp,nPer);
        NTmp=size(tmp);
        tmp = reshape(tmp,[prod(NTmp(1:end-1)),NTmp(end)]);
        
        %         tic;
        tmpH = tmp'*tmp;
        [V,~,~] = svd(tmpH,0);
        tmp = tmp*V;
        tmp = reshape(tmp(:,1:Nc),[NTmp(1:end-1),Nc]);
        %         toc;
        
        tmp = permute(tmp, nPer);
        dataC = dimData(dataC,tmp, slDim,i,'w'); % Write data along 'slDim'
        
    end
    Y = dataC;
end
clear dataC;

end

