function [xHat, RMSE] = GD_SENSE(y, x0, p, method)
%ISTA Summary of this function goes here
%   Detailed explanation goes here
Iteration = p.iteration;
x = x0;
A = p.A; At = p.At;
stepsize = 1/p.L1;
z = x0; t = 1;
RMSE = zeros(1, Iteration);
for i = 1:Iteration
    if method == 1 %GD
        x = x - stepsize*wGradient(y,x,A,At); % gradient decent
    elseif method == 2 % Fast GD
        x_new = z - stepsize*wGradient(y,z,A,At); % gradient decent
        t_new = (1 + sqrt(1 + 4*t^2))/2;
        beta = (t - 1)/t_new;
        z = x_new + beta*(x_new - x);
        x = x_new;
        t = t_new;
    elseif method == 3
        x_new = z - stepsize*wGradient(y,z,A,At); % gradient decent
        t_new = (1 + sqrt(1 + 4*t^2))/2;
        beta = (t - 1)/t_new;
        gamma = t/t_new;
        z = x_new + beta*(x_new - x) + gamma*(x_new - z);
        x = x_new;
        t = t_new;
    end
        RMSE(i) = norm(x(:)-p.xRef(:))/norm(p.xRef(:));
        if mod(i,10) == 0
            disp(['Iteration ' num2str(i) ',' 'Method:' num2str(method), ', RMSE:' num2str(RMSE(i)) ]);
        end
end
xHat = x;

end

