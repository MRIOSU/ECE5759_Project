function display_coil_maps(S)
CH = size(S,4);
figure;
for ch_idx = 1:CH
subplot(2,CH/2,ch_idx);
imagesc(abs(squeeze(S(:,:,1,ch_idx))));
title(['ch' num2str(ch_idx) ' (Magnitude)']);
end

figure;
for ch_idx = 1:CH
subplot(2,CH/2,ch_idx);
imagesc(angle(squeeze(S(:,:,1,ch_idx))));
title(['ch' num2str(ch_idx) ' (Phase)']);
end

end