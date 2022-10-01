% MECH 559 - K. Al Handawi
% McGill University
% Example 4.16 in Papalambros and Wilde

clear
clc
close all
clf

% Visualization of the objective and iterates
x1 = -2:.1:5;
x2=x1;
for i=1:length(x1)
    for j = 1:length(x2)
        f(i,j) = 4*x1(j)^2 + 3*x1(j)*x2(i) + x2(i)^2;
    end
end

% 3D plot
figure(1)
mesh(x1,x2,f)

% contour plot
figure(2)
clf
V = 0:1:20; 
cs = contour(x1,x2,f,V); 
xlim([-2,5])
ylim([-2,5])
hold on

% Choose method (gradient = 1, Newton = 2)
method = 1;

% Ask for intitial guess
xold = [-2 -2]'; 
xold = [5 2]';
xold = [5,5]';

% algorithmic parameters
my_epsilon = .001;
kmax = 100;
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
    g(1) = 8*x(1) + 3*x(2);
    g(2) = 3*x(1) + 2*x(2);
end

function H = hessian(x)
    H = [8 3
         3 2];
end

function f = obj(x)
    f = 4*x(1)^2 + 3*x(1)*x(2) + x(2)^2;
end