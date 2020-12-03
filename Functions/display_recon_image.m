function display_recon_image(xHat_Ref, xHat, select_fr)
     
     xHat_Ref = flip(permute(xHat_Ref, [2 1 3 4 5]),1);
     xHat = flip(permute(xHat, [2 1 3 4 5 6]),1);
     xHat_Ref = padarray(xHat_Ref,[0 1 0 0 0], 0, 'both');
     xHat = padarray(xHat,[0 1 0 0 0 0], 0, 'both');
     xHat_tmp = xHat_Ref;
     error_tmp = zeros(size(xHat_Ref));
     for i = 1:size(xHat,6)
         if i == 1 
             error_scale = 1; 
         else
             error_scale = 10; 
         end
         xHat_tmp = cat(2, xHat_tmp, xHat(:,:,:,:,:,i));
         error_tmp = cat(2, error_tmp, error_scale*(xHat(:,:,:,:,:,i) - xHat_Ref));
     end
     xHat_error = cat(1, xHat_tmp, error_tmp);
     figure; imagesc(abs(squeeze(xHat_error(:,:,:,:,select_fr))), [0, 0.3*max(abs(xHat_Ref(:)))]); axis off; colormap(gray); axis image;
end

