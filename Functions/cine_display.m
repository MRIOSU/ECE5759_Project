function  cine_display(xHat)
%CINE_DISPLAY Summary of this function goes here
%   Detailed explanation goes here
         for fr = 1:size(xHat,5)
             imagesc(abs(xHat(:,:,1,1,fr)), [0,0.3*max(abs(xHat(:)))]); colormap(gray); axis off; title(['frame:' num2str(fr)]); pause(0.07);
         end
end

