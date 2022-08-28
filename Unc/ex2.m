% MECH 559 - M. Kokkolaras
% McGill University
% Example 4.16 in Papalambros and Wilde

clear
clc
%close all
clf

x1 = -2:.1:5;
%x1 = 0:.01:2;
x2=x1;
for i=1:length(x1)
    for j = 1:length(x2)
        f(i,j) = 4*x1(j)^2 + 3*x1(j)*x2(i) + x2(i)^2;
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
xold = [5 2]';
xold = [5,5]';
plot(xold(1),xold(2),'ro')
disp(['objective function value = ',num2str(ex2obj(xold))])

my_epsilon = .001;
kmax = 100;

my_continue = 0;
k = 0;
while my_continue == 0
    k = k + 1
    if method == 1
       % Gradient method
       %pick alpha or use exact line search by commenting accordingly
       alpha = .1;
       %alpha = ex2grad(xold)*ex2grad(xold)'/ ...
               (ex2grad(xold)*ex2hessian(xold)*ex2grad(xold)')
       xnew = xold - alpha*ex2grad(xold)'
    else
       % Newton's method 
       alpha = 1;
       xnew = xold - alpha*inv(ex2hessian(xold))*ex2grad(xold)'
    end
    disp(['norm of gradient = ', num2str(norm(ex2grad(xnew)))])
    disp(['objective function value = ',num2str(ex2obj(xnew))])
    plot(xnew(1),xnew(2),'ro')
    if norm(ex2grad(xnew)) <= my_epsilon
        my_continue = 1;
    end
    xold = xnew;
    if k > kmax
       my_continue = 1;
       disp('maximum number of iteration reached')
   end
end