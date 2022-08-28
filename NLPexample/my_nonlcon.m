function [gineq,geq] = my_nonlcon(x)
%function that evaluates general nonlinear constraints
%returns inequality and equality constraints (negative null formulation)

%initialization
gineq=[];
geq=[];

%evaluate constraints
gineq(1) = x(2) - 1/5*x(1)^2;
 
return

% in case of simulation:
% pass x through an input file
% call/run simulation
% read results 
% gineq = ...
% geq = ...