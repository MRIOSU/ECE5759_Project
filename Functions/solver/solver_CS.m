function [xHat, RMSE] = solver_CS(y, x0, p, method)
%ISTA Summary of this function goes here
%   Detailed explanation goes here
Iteration = p.iteration;
x = x0;
A = p.A; At = p.At;
stepsize = 1/p.L1;
z = x0; t = 1;
theta = 1; gamma = 1; w = x; z = x;
lambda = p.lambda;
RMSE = zeros(1, Iteration);
dim_phase = 5;
for i = 1:Iteration
    if method == 1 % ISTA
        x = x - stepsize*wGradient(y,x,A,At); % gradient decent
        x = shrink(x, lambda/p.L1); % soft thresholding
    elseif method == 2 % FISTA
        x_new = z - stepsize*wGradient(y,z,A,At); % gradient decent
        x_new = shrink(x_new, lambda/p.L1); % soft thresholding
        t_new = (1 + sqrt(1 + 4*t^2))/2;
        z = x_new + (t - 1)/t_new*(x_new - x);
        x = x_new;
        t = t_new;
    elseif method == 3 % POGM
        if i < Iteration - 1
            theta_new = (1 + sqrt(1+4*theta^2))/2;
        else
            theta_new = (1 + sqrt(1+8*theta^2))/2;
        end
        gamma_new = 1/p.L1*(2*theta + theta_new -1)/theta_new;
        w_new = x - stepsize*wGradient(y,x,A,At);
        z_new = w_new + (theta - 1)/theta_new*(w_new - w) + theta/theta_new*(w_new - x) ...
            + (theta - 1)/p.L1/gamma/theta_new*(z - x);
        x = shrink(z_new, lambda*gamma_new);
        w = w_new;
        gamma = gamma_new;
        z = z_new;
        theta = theta_new;
        
    end 
        xHat = ifftc(x, dim_phase);
        RMSE(i) = norm(xHat(:)-p.xRef(:))/norm(p.xRef(:));
        if mod(i,10) == 0
            disp(['Iteration ' num2str(i) ',' 'Method:' num2str(method), ', RMSE:' num2str(RMSE(i)) ]);
        end
end


end

