% MECH 559 - K. Al Handawi
% McGill University
% Example 4.18 in Papalambros and Wilde

clear
clc
close all
clf

% Visualization of the objective and iterates
x = -2:.1:5;
f = x.^4 - 32*x;

figure(1), clf
plot(x,f)
xlim([-2,5])
hold on

% Choose method (gradient = 1, Newton = 2)
method = 2;

% Ask for intitial guess
xold = -3; 

% algorithmic parameters
alpha = .01;
my_epsilon = .1;
kmax = 1000;
iteract = 1; % turn this off to avoid having to hit ENTER on every iteration
verbose = 1; % turn this off to avoid having printing every iteration
plot_progress = 1; % turn this off to avoid plotting every iteration

% optimization
[x_opt,f_opt] = optimize(@obj,@grad,@hessian,method,xold,my_epsilon,kmax,iteract,verbose,plot_progress);

% solution
disp([' ']) % newline
disp(['minimum function value = ',num2str(obj(x_opt))])
disp(['minimizer value = ',sprintf('%0.5f, ',x_opt)])

%% User defined objective, gradient, and hessian
function g = grad(x)
   g = 4*x(1)^3 - 32;
end

function H = hessian(x)
   H = 12*x^2;
end

function f = obj(x)
   f = x^4 -32*x;
end