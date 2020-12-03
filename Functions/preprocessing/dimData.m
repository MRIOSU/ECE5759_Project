function [x] = dimData(x, y, dim, i, a)

% read/write data along the dimension 'dim' at index 'i'

if strcmp(a, 'r')
    if dim==1
        x = x(i,:,:,:,:,:,:);
    elseif dim==2
        x = x(:,i,:,:,:,:,:);
    elseif dim==3
        x = x(:,:,i,:,:,:,:);
    elseif dim==4
        x = x(:,:,:,i,:,:,:);
    elseif dim==5
        x = x(:,:,:,:,i,:,:);
    elseif dim==6
        x = x(:,:,:,:,:,i,:);
    elseif dim==7
        x = x(:,:,:,:,:,:,i);
    else
        error('only works for upto 7 dimensions');
    end

elseif strcmp(a,'w')
    if dim==1
        x(i,:,:,:,:,:,:) = y;
    elseif dim==2
        x(:,i,:,:,:,:,:) = y;
    elseif dim==3
        x(:,:,i,:,:,:,:) = y;
    elseif dim==4
        x(:,:,:,i,:,:,:) = y;
    elseif dim==5
        x(:,:,:,:,i,:,:) = y;
    elseif dim==6
        x(:,:,:,:,:,i,:) = y;
    elseif dim==7
        x(:,:,:,:,:,:,i) = y;
    else
        error('only works for upto 7 dimensions');
    end
    
end