function [est] = k_update(fs,xs,k_old)
% update the Lipschitz constant estimate
alpha = 0.01; k_max = k_old;
    for i=1:length(fs)-1
           k = abs(fs(i)-fs(end))/norm(xs(:,i)-xs(:,end));
           k_max = max(k,k_max);
    end
    i_t = ceil(log(k_max)/log(1+alpha));
    est = (1+alpha)^i_t;
end

