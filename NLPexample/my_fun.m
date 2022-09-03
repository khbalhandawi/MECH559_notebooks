function f = my_fun(x)
% function to evaluate objective function

f=1e6;
f = x(1)^2 + (x(2)-1)^2;
%f=-f;

return

% in case of simulation:
% pass x through an input file
% call/run simulation
% read results 
% f = ...