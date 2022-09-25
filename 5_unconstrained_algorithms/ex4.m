% MECH 559 - M. Kokkolaras
% McGill University
% Example 4.19 in Papalambros and Wilde


clear
clc
clf

x1 = -5:.1:5;
x2=x1;
for i=1:length(x1)
    for j = 1:length(x2)
        f(i,j) = 1/3*x1(j)^3 + x1(j)*x2(i) + 0.5*x2(i)^2 + 2*x2(i) - 2/3;
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
plot(2,-4,'r+'), plot(-1,-1,'r+')

% Choose method (gradient = 1, Newton = 2)
method = 1;

% Ask for intitial guess
%xold = input('Type initial guess as a column vector');
xold = [1 1]'; 
%xold = [-3 2]';
plot(xold(1),xold(2),'ro')
disp(['objective function value = ',num2str(ex4obj(xold))])
%pause
my_epsilon = .001;
kmax = 100;

my_continue = 0;
k = 0;
while my_continue == 0
    k = k + 1
    if method == 1
       % Gradient method
       alpha = ex4grad(xold)*ex4grad(xold)'/ ...
               (ex4grad(xold)*ex4hessian(xold)*ex4grad(xold)')
       xnew = xold - alpha*ex4grad(xold)'
    else
       % Newton's method 
       alpha = 1;
       xnew = xold - alpha*inv(ex4hessian(xold))*ex4grad(xold)'
    end
    disp(['norm of gradient = ', num2str(norm(ex4grad(xnew)))])
    disp(['objective function value = ',num2str(ex4obj(xnew))])
    plot(xnew(1),xnew(2),'ro')
    pause
    if norm(ex4grad(xnew)) <= my_epsilon
        my_continue = 1;
    end
    xold = xnew;
    if k > kmax
       my_continue = 1;
       disp('maximum number of iteration reached')
   end
end