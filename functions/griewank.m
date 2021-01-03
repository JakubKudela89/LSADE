function [val] = griedwank(x)
D = length(x); 
val = 1+(1/4000)*sum(x.^2)-prod(cos(x./sqrt(1:D)));
end