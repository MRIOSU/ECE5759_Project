function Y = noise_prewhiting(data, param)
%input: data [RO E1 E2 CHA PHA SET SLC]
%       param (optinal)
%          param.BW   bandwidth for k-space
%          param.noise [RO, CHA, num_lines] 

N = size(data);
if (nargin > 1 && ~isempty(param.noise))
    nScan = param.noise;
    nScan = squeeze(nScan).* sqrt(param.noise_dwelltime_us/param.acquisition_dwelltime_us);
    nScan = permute(nScan, [2,1,3]);
    nScan = reshape(nScan, [size(nScan,1), numel(nScan)/size(nScan,1)]);
    tmp = 1/(size(nScan,2)-1) * (nScan*nScan'); % covriance matrix
    tmp = (tmp + tmp')/2; % To make the matrix complex hermitian
    L = chol(tmp, 'lower');
    data =  permute(data, [4,1:3,5:numel(N)]); % [CH, RO, E1, E2, PHS, SET, SLC]
    data = reshape(data, [size(data,1), numel(data)/size(data,1)]); % [CH, *]
    data = L\data;
    data = reshape(data, N([4,1:3,5:numel(N)])); % [CH, RO, E1, E2, PHS, SET, SLC]
    Y = permute(data, [2:4,1,5:numel(N)]); % [RO, E1, E2, CH, PHS, SET, SLC]
    
    % check the noise estimated from the k-space after noise-prewhiting
    tmpStd = Y([1:min(16,round(size(Y,1)/8)), end-min(16,round(size(Y,1)/8))+1:end],:,:,:,:,:,:,:,:,:);
    tmpStd = tmpStd(:, [1:round(size(Y,2)/5), end-round(size(Y,2)/5)+1:end],:,:,:,:,:,:,:,:);
    tmpStd = tmpStd(:, :, [1:max(1,round(size(Y,3)/5)), end-round(size(Y,3)/5)+1:end],:,:,:,:,:,:,:,:);
    mStd = zeros(size(Y,4),1);    
    for i = 1:size(Y,4)
        tmp = tmpStd(:,:,:,i,:,:);  tmp = std(tmp(tmp~=0));
        mStd(i) = tmp; % Noise in the ith coil
    end
    disp('Noise estimate from the k-space after noise-prewhiting:')
    disp(num2str(mStd'));
    
else
    warning('Noise is estimated from the k-space data');
    tmpStd = data([1:min(16,round(size(data,1)/8)), end-min(16,round(size(data,1)/8))+1:end],:,:,:,:,:,:,:,:,:);
    tmpStd = tmpStd(:, [1:round(size(data,2)/5), end-round(size(data,2)/5)+1:end],:,:,:,:,:,:,:,:);
    tmpStd = tmpStd(:, :, [1:max(1,round(size(data,3)/5)), end-round(size(data,3)/5)+1:end],:,:,:,:,:,:,:,:);
    mStd = zeros(size(data,4),1);
    Y = zeros(size(data));
    for i = 1:size(data,4)
        tmp = tmpStd(:,:,:,i,:,:);
        tmp = std(tmp(tmp~=0));
        mStd(i) = tmp; % Noise in the ith coil
        Y(:,:,:,i,:,:,:) = data(:,:,:,i,:,:,:)/tmp;
    end
    disp(num2str(mStd'));
    % tmp = prctile(sort(mStd), 25)/1.2; % pick a coil close to the bottom in terms of signal power
    % Y = data*nStd/tmp;
end

end

