function [val] = rosenbrock(x)
D =length(x); val = 0;
for i=1:D-1
    val = val + 100*(x(i+1)-x(i)^2)^2 + (1-x(i))^2;
end
end