function [val] = fub(f_xi,xi,k,x)
% Lipschitz-based surrogate value at point x
t = length(f_xi);
vals = zeros(size(f_xi));
for i=1:length(f_xi)
    vals(i) = min(f_xi(i) + k*norm(xi(:,i)-x));
end
val = min(vals);
end