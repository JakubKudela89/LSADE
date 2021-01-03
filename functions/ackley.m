function [val] = ackley(x)
D = length(x); 
val = -20*exp(-0.2*sqrt((1/D) * sum(x.^2))) - exp((1/D) * sum(cos(2*pi*x))) + 20 + exp(1);
end