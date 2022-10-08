clear all
clc
close all
format short
format compact

% % corn/oats example
% A = [1 1 1 0; 2 1 0 1];
% b = [ 240 320 ]';
% c = [ -40 -30 0 0];

% % kauserwise's example
% A = [10 20 1 0; 8 8 0 1];
% b = [ 120 80 ]';
% c = [ -12 -16 0 0];

% Kokkolaras' example (with negative values)
A = [2 3 -2 -3 1 0 0; -5 -2 5 2 0 1 0; -2 7 2 -7 0 0 1];
b = [10; -2; 8];
c = [1 -1 -1 1 0 0 0];

% Using MATLAB's linprog
LB = zeros(1,size(A,2));
UB = Inf*ones(1,size(A,2));
fprintf("MATLABs linprog solution\n")
xopt = linprog(c,[],[],A,b,LB,UB)

stop = false;
k = 0;
cf = c;
fprintf("\n==================" + ...
    "\nSimplex iterations\n")
while ~stop
    [x_opt,stop,A,b,c,k] = simplex(A,b,c,k,true);
    f_opt = cf*x_opt
    fprintf("\n==================\n")
    pause
end