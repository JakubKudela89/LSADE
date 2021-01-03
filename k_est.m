function [est] = k_est(fs,xs)
% estimate the Lipschitz constant
alpha = 0.01; k_max = 0;
    for i=1:length(fs)-1
       for j=i+1:length(fs)
           k = abs(fs(i)-fs(j))/norm(xs(:,i)-xs(:,j));
           k_max = max(k,k_max);
       end
    end
    i_t = ceil(log(k_max)/log(1+alpha));
    est = (1+alpha)^i_t;
end