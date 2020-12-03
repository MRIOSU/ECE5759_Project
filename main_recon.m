% This is the Matlab script for ECE 5759 final project
% Last modified: 12-03-2020 by Chong Chen (Chong.Chen@osumc.edu)

close all;
restoredefaultpath
addpath(genpath('./Functions'));

% parmaters for coil sensitivity map
p.mthd   = 1; %'1' espirit, '2' Walsh
p.fil    = [6,6,6];  % size of kernal for eSPIRiT or size of filter for 'Walsh'; use [6,6,6] for eSPIRiT and [3,1,1] for Walsh
p.ACSsz  = [128,96,30]; % size of the k-space block used for eSPIRiT
p.eSRT   = 0.95;
p.fixeSRTsig = []; % fixed sigular value threshold (just for debug, set empty when running)
p.eSRTsig = 0.5; % sigular value threshold setted manually
p.eSmaps = 1; %number of Espirit sensitivity maps
p.avgPhs = 1; % Assign the phase to time-average image. 1: yes, 0: no

[fileName,dirName,FilterIndex] = uigetfile('*.h5','MultiSelect', 'on');
for k = 1:size(cellstr(fileName),2)
    if iscell(fileName)
        filename = fileName{k};
    else
        filename = fileName;
    end
    
    %% read the data
    [kData, param] = read_ocmr([dirName filename]); % Kdata is orgnazide as {'kx'  'ky'  'kz'  'coil'  'phase'  'set'  'slice'  'rep'  'avg'}
    kData_org = kData; param_org = param; % backup the data
    
    %% preprocessing
    % remove oversampling
    kData = remove_RO_oversamling(kData);
    param.noise = remove_RO_oversamling(param.noise);
    % noise prewhiting
    kData = noise_prewhiting(kData, param);
    % coil compression
    nCH = 12;
    kData = coil_compression(kData, nCH);
    
    %% define the operator A
    % calculate coil sensitivity maps
    tic;[S,x0]  = coilSen(kData, p);toc
    PHS = size(kData, 5);
    x0  = repmat(x0,[1 1 1 1 PHS]);%Insert frame dim
    x0  = x0 + 0.010*max(abs(x0(:)))*randn(size(x0)); % inital image guess
    [y, p.A, p.At]= sensor_operator(kData, S);
    
    %% recon fully sampled dataset A^T(y) (reference)
    xHat_IFFT = p.At(y);
    figure;
    for rep = 1:3
        cine_display(xHat_IFFT);
    end
    p.xRef = xHat_IFFT;
    
    %% downsampling
    E1 = size(kData,2); FR = size(kData,5);
    R = 4; %retrospective downsampling factor
    
    %% uniform sampling (parall imaging)
    for R = 2:2:6
        uniform_sampling = zeros([1,E1,1,1,FR]);
        for fr = 1:size(uniform_sampling,5)
            uniform_sampling(1,mod(fr,R)+1:R:end,:,:,fr) = 1;
        end
        figure;subplot(1,3,1:2);imagesc(repmat(rot90(squeeze(uniform_sampling(:,:,1,1,1,1,1))),[1 size(kData,1)]));  colormap(gray);xlabel('kx','FontSize',14); ylabel('ky','FontSize',14)
        subplot(1,3,3);imagesc(squeeze(uniform_sampling));  colormap(gray); xlabel('t','FontSize',14); ylabel('ky','FontSize',14)
        title('Unifrom sampling pattern for paramel imaging (R = 4)')
        DataIn_uniform = bsxfun(@times, kData, uniform_sampling);
        if R == 2
            p.iteration = 30;
        else
            p.iteration = 150;
        end
        % SENSE reconstruction
        [y, p.A, p.At]= sensor_operator(DataIn_uniform, S);
        L1 = powerIter(p.A, p.At, size(x0));
        p.L1 = L1*2.05;
        RMSE_uniform = zeros(p.iteration, 3);
        xHat_uniform = zeros([size(x0),3]);
        for method = 1:3
            [xHat_uniform(:,:,:,:,:,method), RMSE_uniform(:,method)] = GD_SENSE(y, randn(size(x0)), p, method);
        end
        iter_step = 3; iter_x = 1:iter_step:size(RMSE_uniform,1);
        figure; plot(iter_x, RMSE_uniform(1:iter_step:end,1),'*-'); hold on; plot(iter_x, RMSE_uniform(1:iter_step:end,2),'d-'); plot(iter_x, RMSE_uniform(1:iter_step:end,3),'o-'); lgd = legend('GM','FGM','OGM');
        lgd.FontSize = 14; xlabel('Iteration#', 'FontSize', 14); ylabel('NRMSE', 'FontSize', 14); title(['R = ' num2str(R)], 'FontSize', 14);
        grid on; set(gcf, 'Position', [680   698   367   280]);
        figure;
        for rep = 1:2
            cine_display(xHat_uniform(:,:,:,:,:,3));
        end
        xHat_inv = R*p.At(y);
        display_recon_image(xHat_IFFT, cat(6,xHat_inv, xHat_uniform), 11);
    end
    
    
    
    %% random sampling (parallal imaging + compressive sensing)
    for R = 4:2:8
        sampling_pattern = GoFIX(E1,FR,R); %[E1, FR, SET]
        random_sampling = logical(permute(sampling_pattern,[4 1 5 6 2 3]));
        DataIn_random = bsxfun(@times, kData, random_sampling);
        figure;subplot(1,3,1:2);imagesc(repmat(rot90(squeeze(random_sampling(:,:,1,1,1,1,1))),[1 size(kData,1)]));  colormap(gray);xlabel('kx','FontSize',14); ylabel('ky','FontSize',14)
        subplot(1,3,3);imagesc(squeeze(random_sampling));  colormap(gray); xlabel('t','FontSize',14); ylabel('ky','FontSize',14)
        
        % define sparsity transformation
        dim_phase = 5;
        [y, p.A, p.At]= sensor_operator(DataIn_random, S);
        p.A = @(x) p.A(ifftc(x,dim_phase));
        p.At = @(y) fftc(p.At(y),dim_phase);
        L1 = powerIter(p.A, p.At, size(x0));
        p.L1 = L1*2.05;
        p.iteration = 150;
        p.lambda = 0.5*10^5*scale_factor;
        RMSE_random = zeros(p.iteration, 3);
        xHat_random = zeros([size(x0),3]);
        for method = 1:3
            tic;[xHat_random(:,:,:,:,:,method), RMSE_random(:,method)] = GD_CS(y, randn(size(x0)), p, method);toc;
            disp('-------------')
        end
        iter_step = 3; iter_x = 1:iter_step:size(RMSE_random,1);
        figure; plot(iter_x, RMSE_random(1:iter_step:end,1),'*-'); hold on; plot(iter_x, RMSE_random(1:iter_step:end,2),'d-'); plot(iter_x, RMSE_random(1:iter_step:end,3),'o-'); lgd = legend('ISTA/PGM','FISTA/PFGM','POGM');
        lgd.FontSize = 14; xlabel('Iteration#', 'FontSize', 14); ylabel('NRMSE', 'FontSize', 14); title(['R = ' num2str(R)], 'FontSize', 14);
        grid on; set(gcf, 'Position', [680   698   367   280]);
        figure;
        for rep = 1:2
            cine_display(xHat_random(:,:,:,:,:,3));
        end
        xHat_inv = p.At(y);
        display_recon_image(xHat_IFFT, cat(6,xHat_inv, xHat_random), 11);
    end
    
    % find optimal lambda
%     p.iteration = 150;
%     RMSE_CS = zeros(p.iteration, 9);
%     i = 1;
%     for scale_factor = 2.^[-4:4]
%         p.lambda = 0.5*10^5*scale_factor;
%         % ISTA reconstruction
%         method = 1;
%         [xHat_ISTA, RMSE_CS(:,i)] = GD_CS(y, randn(size(x0)), p, method);
%         i = i + 1;
%     end
    
end
