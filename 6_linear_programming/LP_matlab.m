clear all
clc
close all
format short
format compact

%% Using MATLAB's linprog
c =   [-1, 1];
%      ─┬  ─┬
%       │   └┤ Coefficient for x1
%       └────┤ Coefficient for x2
A_ineq =   [ 2,  3;...  % constraint g1 LHS
            -5, -2;...  % constraint g2 LHS
            -2,  7;...  % constraint g3 LHS
             1,  0;...  % constraint g4 LHS
             0, -1];    % constraint g5 LHS
b_ineq =   [10,...  % constraint g1 RHS
            -2,...  % constraint g2 RHS
             8,...  % constraint g3 RHS
             5,...  % constraint g4 RHS
             5];    % constraint g5 RHS

LB = [-Inf, -Inf];
%      ─┬     ─┬
%       │      └┤ lower bound for x1
%       └────┤ lower bound for x2
UB = [Inf, Inf];
%     ─┬   ─┬
%      │    └┤ upper bound for x1
%      └────┤ upper bound for x2

fprintf("MATLABs linprog solution\n")
[x_opt,f_opt] = linprog(c,A_ineq,b_ineq,[],[],LB,UB)

%% Using simplex.m
% Kokkolaras' example (with negative values)
A = [...
    2 3 -2 -3 1 0 0 0 0; ...
    -5 -2 5 2 0 1 0 0 0; ...
    -2 7 2 -7 0 0 1 0 0; ...
    0 1 0 -1 0 0 0 -1 0; ...
    1 0 -1 0 0 0 0 0 1];
b = [10 -2 8 -5 5]'; % remember to make b a column vector
c = [-1 1 1 -1 0 0 0 0 0]; % this is c transpose

stop = false;
k = 0;
cf = c;
fprintf("\n==================" + ...
    "\nSimplex iterations\n")
while ~stop
    input("hit ENTER to continue")
    [x_opt,stop,A,b,c,k] = simplex(A,b,c,k,true);
    z_opt = x_opt(1:2) - x_opt(3:4)
    f_opt = cf*x_opt
    fprintf("\n==================\n")
end

%% Other examples

% % corn/oats example
% A = [1 1 1 0; 2 1 0 1];
% b = [ 240 320 ]';
% c = [ -40 -30 0 0];

% % kauserwise's example
% A = [10 20 1 0; 8 8 0 1];
% b = [ 120 80 ]';
% c = [ -12 -16 0 0];