function xHat = ISTA(y, x0, p)
%ISTA Summary of this function goes here
%   Detailed explanation goes here

Iteration = p.iteration;
x = x0;
A = p.A; At = p.At;
stepsize = 1/p.L1;
lambda = p.lambda;
for i = 1:Iteration
    x = x - stepsize*wGradient(y,x,A,At); % gradient decent
%     x = shrink(x, lambda/p.L1); % soft thresholding
    if mod(i,10) == 0
        disp(['Iteration ' num2str(i) ',' ]);
    end
end
% dim_Phase = 5;
% xHat = fftc(x, dim_Phase);
xHat = x;

end

