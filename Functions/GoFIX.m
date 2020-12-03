function [samp] = GoFIX(PE, FR, R)


%% Essential paramters
n = round(PE/R);   % Number of PE lines per frame
% FR = 48;   % Frames
% PE  = 144;   % Size of of PE
% special cart parameters
k  = 3;   % k>=1; % k=1 uniform; k>1 variable density profile; larger k means flatter top (default: 3, range: 1-10, precision: 0.1)
s  = 0.9;   % s>=0; % largers s means more higher sampling density in the middle (default: 1, range: 0-10, precision: 0.1)
PF = 0;     % for partial fourier; discards PF samples from one side (default: 0, range: 0-floor(n/2), precision: 1);


%% Other parameters
ir = 1;     % ir = 1 or 2 for golden angle, ir > 2 for tiny golden angles
gr = (1+sqrt(5))/2; % golden ratio F_{n+1}/F_{n}
ga = 1/(gr+ir-1); %(1-1/gr); % golden angle, sqrt(2) works equally well
N  = n*FR; % Total number of k-space line collected
R  = PE/(n+PF); % Acceleration
s  = max(1, s*(R).^(1/3));     % % To keep the default value of s close to 1
PES = ceil(PE * 1/s); % Size of shrunk PE grid


%% Size of the smaller pseudo-grid which after stretching gives the true grid size
vd = (PE/2-PES/2)/((PES/2)^k); %(1/2*(n-PES)/max(tmp)); % location specific displacement

%% Let's populate the grid;
samp  = zeros(PE, FR); % sampling on PE-t grid

figure;
kk=1;
PEInd = zeros((n-PF)*FR,1); % Final PE index used for MRI
v0 = (1/2+1e-10:PES/(n+PF):PES+1/2-1e-10); % Start with uniform sampling for each frame

% E = 2; % These few line are for PC-MRI with E>1
% for e=0:E-1
% v0 = v0+e*PES/(E*(n+PF)); % Start with uniform sampling for each frame
% kk=e+1;
for j=1:FR
    v = rem((v0 + (j-1)*PES/(n+PF)*ga)-1, PES) + 1; % In each frame, shift by golden shift of PES/TR*ga
    v = v - PES.*(v>=(PES+0.5));
    
    if rem(PE,2)==0 % if even, shift by 1/2 pixel
        vC = v - vd*sign((PES/2+1/2)-v).*(abs((PES/2+1/2)-v)).^k + (PE-PES)/2 + 1/2;%(ctrn - ctrnS);
        vC = vC - PE.*(vC>=(PE+0.5));
    elseif rem(PE,2)==1 % if odd don't shift
        vC = v - vd*sign((PES/2+1/2)-v).*(abs((PES/2+1/2)-v)).^k + (PE-PES)/2;%(ctrn - ctrnS);
    end
    vC = round(sort(vC));
    vC(1:PF) = [];
    
    if rem(j,2) == 1
        PEInd((j-1)*n + 1 : j*n) = vC;
    else
        PEInd((j-1)*n + 1 : j*n) = flip(vC);
    end
    
    samp(vC, j) = samp(vC, j)+ kk;
    imagesc(samp); xlabel('frames'); ylabel('PE'); %axis('image'); %colormap(c); %colorbar;
    pause(1e-3);
end
% end

%%
% PEInd = find(samp~=0);
% figure; imagesc(abs(fftshift(fft(ifftshift(samp))))); axis('image'); title('PSF');
% figure; imagesc(abs(fftshift(fft2(ifftshift(samp))))); axis('image'); title('PSF');
% 
% figure; plot(sum(samp,2)); xlabel('PE index'); ylabel('freq'); title('Time-averaged distribution');
% % figure; imagesc(samp); colormap(hot);%colorbar; %axis('image');
% figure; plot(PEInd(1:end),'*-'); xlabel('line index'); ylabel('PE index'); % PE order
% figure; hist(diff(PEInd),32); title('Distribution of k-space jumps')







