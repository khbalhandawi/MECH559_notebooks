function [f,dfdx] = my_funED(x)
% function to evaluate objective function

f = x(1)^2 + (x(2)-1)^2;
dfdx(1) = 2*x(1);
dfdx(2) = 2*(x(2)-1);

% in case of simulation:
% write input file as function of x
% call and run simulation
% read results as output
%f = some function of the output