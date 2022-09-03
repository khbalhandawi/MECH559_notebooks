function [gineq,geq,dgineqdx,dgeqdx] = my_nonlconED(x)
%function that evaluates general nonlinear constraints
%returns inequality and equality constraints (negative null formulation)

%initialization
gineq=[];
geq=[];
dgineqdx=[];
dgeqdx=[];

%evaluate constraints and derivatives
gineq(1) = x(2) - 1/5*x(1)^2;
dgineqdx = [-2/5*x(1); 1];
%dgineqdx(2) = 1;


