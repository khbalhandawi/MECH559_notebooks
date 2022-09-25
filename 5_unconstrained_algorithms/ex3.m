% MECH 559 - M. Kokkolaras
% McGill University
% Example 4.18 in Papalambros and Wilde

clear
clc
close all
clf

x = -2:.1:5;
f = x.^4 - 32*x;

figure(1), clf
plot(x,f), hold on

% Choose method (gradient = 1, Newton = 2)
method = 2;

% Ask for intitial guess
%xold = input('Type initial guess as a column vector');
xold = -3; 
plot(xold,ex3obj(xold),'ro')
disp(['objective function value = ',num2str(ex3obj(xold))])

alpha = .01;
my_epsilon = .1;
kmax = 1000;

my_continue = 0;
k = 0;
while my_continue == 0
    k = k + 1
    if method == 1
       % Gradient method
       alpha = ex3grad(xold)*ex3grad(xold)'/ ...
               (ex3grad(xold)*ex3hessian(xold)*ex3grad(xold)')
       xnew = xold - alpha*ex3grad(xold)'
    else
       % Newton's method 
       alpha = .01;
       xnew = xold - alpha*inv(ex3hessian(xold))*ex3grad(xold)'
    end
    disp(['norm of gradient = ', num2str(norm(ex3grad(xnew)))])
    disp(['objective function value = ',num2str(ex3obj(xnew))])
    plot(xnew,ex3obj(xnew),'ro')
    pause
    if norm(ex3grad(xnew)) <= my_epsilon
        my_continue = 1;
    end
    xold = xnew;
    if k > kmax
       my_continue = 1;
       disp('maximum number of iteration reached')
   end
end