function [f,c] = myexample(x)

f=1e6;
c=[];

f = x(1)^2 + (x(2)-1)^2;
c(1) = x(2) - 1/5*x(1)^2;