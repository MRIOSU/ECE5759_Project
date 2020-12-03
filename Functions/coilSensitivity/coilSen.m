%%
%Input:
%kdata: data in k-space
%       [kx, ky, kz, coil, phase, set, slc] %slc=1
%samp: sampling pattern
%      [kx, ky, kz, coil, phase, set, slc] %slc=1
%      type: boolean, true indicating the point was acquired
%param:
%     .mthd  : method
%     .fil   : order of smoothing matrix, default 9
%     .opt   : type of phase correction

%Output:
%S: sensitivity map
%   [RO, E1, E2, CHA, num_of_coil_maps]
%cI: [RO, E1, E2, num_of_coil_maps]

%%
function [S, cI] = coilSen(kdata, param)

mthd = param.mthd;
fft3Scal = 1/sqrt(size(kdata,1)*size(kdata,2)*size(kdata,3));

    % average over time for each coil
    dims_phase = 5;
    dims_set = 6;
    samp = logical(abs(kdata));
    samp = samp(:,:,:,1,:,:,1); %set dim_CHA, dim_SLC = 1    
    samp_avg = sum(sum(samp,dims_phase),dims_set);    
    kdata_avg = sum(sum(kdata,dims_phase),dims_set)./(repmat(samp_avg,[1,1,1,size(kdata,4)])+eps);
    samp_avg = logical(samp_avg);

if mthd == 1 % espirit
    disp(['Espirit,' num2str(param.eSmaps) ' sensitivity maps']);
%     CalibSize = getCalibSize3D(logical(sum(samp,4)));
    CalibSize = size(samp_avg);
    if numel(CalibSize) == 2
        CalibSize(3) = 1;
    end
    CalibSize = min(CalibSize,param.ACSsz); % specify the calibration region
    [cI,S] = espirit_sens3D(kdata_avg,samp_avg,CalibSize,param);
    [cI,S] = sensCorrect3D(cI, S, param.avgPhs); % Remove phase from time average    
%     S = squeeze(S); cI = squeeze(cI); % pesudo-3D to 2D

elseif mthd == 2 % Walsh
    disp('Walh,1 sensitivity map');    
    if size(kdata_avg,3) == 1 % pesudo-3D to 2D
        fft2Scal = fft3Scal;
        I = kdata2im(squeeze(kdata_avg),squeeze(samp_avg), fft2Scal, param.ACSco); % Find low-res images
        [cI,S] = WalshCoilCombine(I, param);
        [cI,S] = sensCorrect(cI, S, param.avgPhs); % Remove phase from time average 
        S = permute(S,[1 2 4 3]); %RO E1 E2 CH
        cI = permute(cI,[1 2 3]); %RO E1 E2
    else % 3D
        I = kdata2im3D(kdata_avg,samp_avg, fft3Scal); % Find low-res images
        [S,cI] = WalshCoilCombine3D(I, param);
%     [cI,S] = sensCorrect3D_m(cI, S, param.avgPhs); % Remove phase from time average(already in WalshCoilCombine3D)         
    end
    
else
    error('Undefined method for coil sensisitive estimation');
end




% Find low-res image from the data 'd'
function img = kdata2im(d, s, fft2Scal, co)
[CalibSize, ~] = getCalibSize(s);
fil = g2d(size(d(:,:,1)), CalibSize, co);
d = d.* repmat(fil, [1,1,size(d,3)]);
img = 1/fft2Scal * fftshift(fftshift(ifft2(ifftshift(ifftshift(d,1),2)),2),1);



function F = g2d(S1, S2, co)
[x,y] = ndgrid(1:S1(1), 1:S1(2));
cntr=floor(S1/2)+1;
sigx = S2(1)*co(1);
sigy = S2(2)*co(2);
F = 1/(2*pi*sigx*sigy)*exp(-((x-cntr(1)).^2)./(2*sigx^2) - ((y-cntr(2)).^2)./(2*sigy^2));
F = F/max(F(:));




% Find low-res image from the data 'd'
function img = kdata2im3D(d, s, fft3Scal)
% slicer(s,1,1,' ');
sz = size(s); % Size of s

for j = 1 : numel(size(s)) % Go to s, dimension by dimension
    if j == 1
        tmp = squeeze(s(:, floor(end/2)+1, floor(end/2)+1))';
    elseif j == 2
        tmp = squeeze(s(floor(end/2)+1, :, floor(end/2)+1))';
    elseif j == 3
        tmp = squeeze(s(floor(end/2)+1, floor(end/2)+1, :))';
    end
    dis = abs(-sz(j)/2+0.5 : 1 : sz(j)/2-0.5)';
    tmpZ = find(~tmp);
    if ~isempty(tmpZ)
        [mnVal, ~]= min(dis(tmpZ));
        dis(dis>(mnVal-0.5))= 0;
        dis = logical(dis);
    else
        dis = ones(sz(j),1);
    end
    CSize(j) = sum(dis);
end
% crpScl(1) = param.crpScl(1);
% crpScl(2) = param.crpScl(2);
% size(s),
% slicer(s,1,1,'');
% cntr = floor(size(s)/2)+1;
% [indx, indy, indz] = ind2sub(size(s), find(s+1));
% indx = indx - cntr(1);
% indy = indy - cntr(2);
% indz = indz - cntr(3);
% dis = sqrt(crpScl(1)*indx.^2 + indy.^2 + crpScl(2)*indz.^2);
% disTmp = dis(~s);
% [mnVal,~] = min(disTmp);
% dis(dis>mnVal) = 0;
% dis = logical(reshape(dis, size(s)));
% [indx, indy, indz] = ind2sub(size(dis), find(dis));
% indx = 2*max(abs(indx - cntr(1))),
% indy = 2*max(abs(indy - cntr(2))),
% indz = 2*max(abs(indz - cntr(3))),
% slicer(dis,1,1,'');

% CalibSize = getCalibSize3D(s); %% usinng the new calibsize
% CSize = CalibSize;

fil = g3d(size(d(:,:,:,1)), CSize); % Lowpass filter
blk = zeros(size(d(:,:,:,1))); % Central continuous acquired block
blk(1:CSize(1), 1:CSize(2), 1:CSize(3)) = 1;
dShft = floor((size(blk) - CSize)/2);
blk = circshift(blk, [dShft(1), dShft(2), dShft(3)]);
slicer(blk.*fil,1,1,'');

% zSizeS = floor((size(d(:,:,:,1)) - CSize)/2);
% zSizeE = ceil((size(d(:,:,:,1)) - CSize)/2)-1;
% fil(1:zSizeS(1),:,:) = 0;
% fil(end-zSizeE(1):end,:,:) = 0;
% 
% fil(:,1:zSizeS(2),:)       = 0;
% fil(:,end-zSizeE(2):end,:) = 0;
% 
% fil(:,:,1:zSizeS(3))       = 0;
% fil(:,:,end-zSizeE(3):end) = 0;

% slicer(abs(fil),1,1.0, 'calibration window');

d = d .* repmat(blk, [1,1,1,size(d,4)]) .* repmat(fil, [1,1,1,size(d,4)]);
img = zeros(size(d));
for c = 1:size(d,4)
    img(:,:,:,c) = 1/fft3Scal * fftshift(fftshift(fftshift(ifftn(ifftshift(ifftshift(ifftshift(d(:,:,:,c),1),2),3)),3),2),1);
%     slicer(abs(img(:,:,8:10,c)),1,1,'');
end



function F = g3d(S1, S2)
[x,y,z] = ndgrid(1:S1(1), 1:S1(2), 1:S1(3));
cntr=floor(S1/2)+1; % center
sigx = S2(1)/3; % sigma x
sigy = S2(2)/3;%max(S2(2)/1, sigx/3); % sigma y
sigz = S2(3)/3;%max(S2(3)/1, sigx/3); % sigma z
F = 1/(2*pi*sigx*sigy*sigz)*exp(-((x-cntr(1)).^2)./(2*sigx^2) - ((y-cntr(2)).^2)./(2*sigy^2) - ((z-cntr(3)).^2)./(2*sigz^2));
F = F/max(F(:));
% slicer(F,1,1,'');
% figure; imagesc(F); axis('image');

