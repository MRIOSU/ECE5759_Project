function [cI, S] = WalshCoilCombine(I, param)
%
% ------------------------ input ------------------------
% I: coil-by-coil complex images [FE, PE, Coils]. 
% param: input parameters
%        param.fil = filter size used to smoothen E[X'X] matrix
%        param.opt = type of phase correction
       

% ------------------------ output ------------------------
% cI: coil-combined image
% S: Sensitivity maps

% ---------------------- References ----------------------   
% Polarimetric techniques for enhancing SAR imagery, SPIE, 1630: 141-173 1992
% Adaptive Reconstruction of Phased Array MR Imagery, MRM 43:682-690 2000


% Check input variables
if nargin~=1 && nargin~=2
    error('Incorrect number of input arguments');
    
elseif nargin==1 
    param.fil = [3,1];
    param.opt = 2;
end

if ~isfield(param, 'opt')
    param.opt = 2;
elseif param.opt<1 || param.opt>3
    param.opt = 2;
    warning('param.opt is set to 2');
end

if ~isfield(param, 'fil')
    param.fil = [3,1];
elseif param.fil(1)<1
    param.fil = [3,1];
    warning('param.fil is set to 9');
end
    

N = size(I); % Image size
FE = N(1); % Size of frequency encoding direction
PE = N(2); % Size of phase encoding direction
if numel(N) ==2
    Nc = 1;
else
    Nc = N(3); % Number of coils
end
param.fil = 2*(floor(param.fil/2))+1; % To ensure it is odd

cI = zeros(FE,PE);
S  = zeros(FE,PE,Nc);


%% Bulid "covariance" matrix covI
covI = zeros(FE,PE,Nc,Nc);
for k=1:Nc
    covI(:,:,k,k)=I(:,:,k).*conj(I(:,:,k));
    for j=1:k-1
		covI(:,:,k,j)=I(:,:,k).*conj(I(:,:,j));
		covI(:,:,j,k)=conj(covI(:,:,k,j));
    end
end

% Spatially smoothen covI to generate covIs
fil2 = ones(param.fil);%, param.fil);
covIs = zeros(FE,PE,Nc,Nc);
for i = 1:Nc
    for j=1:Nc
%         covIs(:,:,i,j) = conv2(covI(:,:,i,j), fil2, 'same');
            covIs(:,:,i,j) = covI(:,:,i,j);
%         covIs(:,:,i,j) = medfilt2(real(covI(:,:,i,j)), [13,13]) + ...
%                       1j* medfilt2(imag(covI(:,:,i,j)), [13,13]);
%         figure; imagesc(squeeze(imag(covIs(:,:,i,j))));
%         covIs(:,:,i,j) = conv2(real(covI(:,:,i,j)), fil, 'same') + 1j*conv2(imag(covI(:,:,i,j)), fil, 'same');
    end
end

% Note: This way of smoothing covI has been suggested by Peter Kellman and
% reportedly works well.


%% Performing coil combine and senstitivity estimation per Walsh et al.
for i = 1:PE
    for j = 1:FE               
        covIs_2D = reshape(covIs(j, i, :, :), [Nc Nc]);
        [U,~] = eig((covIs_2D + covIs_2D')/2); 
        cI(j,i) = reshape(I(j, i, :), [1 Nc]) * conj(U(:,end));
        S(j,i,:) = reshape(U(:,end), [1 1 Nc]);   
%         S(j,i,:) = reshape(abs(U(:,end)).^1.5 .* (U(:,end)./abs(U(:,end))), [1 1 Nc]);   

        
%         covIs_2D = reshape(covIs(j, i, :, :), [Nc Nc]);
%         [U,D] = eig((covIs_2D + covIs_2D')/2);        
%         loads{1} = U*D;
%         loads{2} = U;
%         [~,loads] = sign_flip(loads,(covIs_2D + covIs_2D')/2);
%         cI(j,i) = reshape(I(j, i, :), [1 Nc]) * conj(loads{2}(:,end));
%         S(j,i,:) = reshape(loads{2}(:,end), [1 1 Nc]);    
    end
end

