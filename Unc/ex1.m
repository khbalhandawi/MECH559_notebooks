% MECH 559 - M. Kokkolaras
% McGill University
% Example 4.15 in Papalambros and Wilde

clear all
clc
%close all
clf

x1 = -5:.1:5;
%x1 = -2:.1:2;
x1 = -10:.1:10;
%x1 = 0:.01:2;
x2=x1;
for i=1:length(x1)
    for j = 1:length(x2)
        f(i,j) = x1(j)^4-2*x1(j)^2*x2(i)+x2(i)^2;
    end
end

figure(1)
mesh(x1,x2,f)

figure(2)
clf
V = 0:1:20; 
cs = contour(x1,x2,f,V); 
clabel(cs)
hold on

% Choose method (gradient = 1, Newton = 2)
method = 2;

% Ask for intitial guess
%xold = input('Type initial guess as a column vector');
xold = [-2 -2]'; 
%xold = [5 2]';
xold = [2,1]';
%xold = [1,0]';
plot(xold(1),xold(2),'ro'), grid on
disp(['objective function value = ',num2str(ex1obj(xold))])

my_epsilon = .000001;
kmax = 1000;

my_continue = 0;
k = 0;
while my_continue == 0
    k = k + 1
    if method == 1
       % Gradient method
       %pick alpha or use exact line search by commenting accordingly
       alpha = 0.01;
       alpha = ex1grad(xold)*ex1grad(xold)'/ ...
               (ex1grad(xold)*ex1hessian(xold)*ex1grad(xold)')
       xnew = xold - alpha*ex1grad(xold)'
    else
       % Newton's method (won't work because of singular Hessian)
       alpha = 1;
       xnew = xold - alpha*inv(ex1hessian(xold))*ex1grad(xold)'
    end
    disp(['norm of gradient = ', num2str(norm(ex1grad(xnew)))])
    disp(['objective function value = ',num2str(ex1obj(xnew))])
    pause
    plot(xnew(1),xnew(2),'ro')
    if norm(ex1grad(xnew)) <= my_epsilon
        my_continue = 1;
    end
    xold = xnew;
    if k > kmax
       my_continue = 1;
       disp('maximum number of iteration reached')
    end
end