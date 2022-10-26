% MECH 559 - M. Kokkolaras
% McGill University
% 

clear all
clc
%close all
clf

x1 = -2:.1:2;
%x1 = -2:.1:2;
x1 = -2:.1:2;
%x1 = 0:.01:2;
x2=x1;
for i=1:length(x1)
    for j = 1:length(x2)
        f(i,j) = x1(j)+(1/x1(j))+x2(i)+(1/x2(i));
    end
end

figure(1)
mesh(x1,x2,f)

figure(2)
clf
V = 0:1:20; 
cs = contour(x1,x2,f,V); 
hold on

% Choose method (gradient = 1, Newton = 2)
method = 1;

% Ask for intitial guess
%xold = input('Type initial guess as a column vector');
xold = [-2 -2]'; 
%xold = [5 2]';
xold = [2,1]';
%xold = [1,0]';
plot(xold(1),xold(2),'ro'), grid on
disp(['objective function value = ',num2str(obj(xold))])

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
       alpha = grad(xold)*grad(xold)'/ ...
               (grad(xold)*hessian(xold)*grad(xold)')
       xnew = xold - alpha*grad(xold)'
    else
       % Newton's method (won't work because of singular Hessian)
       alpha = 1;
       xnew = xold - alpha*inv(hessian(xold))*grad(xold)'
    end
    disp(['norm of gradient = ', num2str(norm(grad(xnew)))])
    disp(['objective function value = ',num2str(obj(xnew))])
    pause
    plot(xnew(1),xnew(2),'ro')
    if norm(grad(xnew)) <= my_epsilon
        my_continue = 1;
    end
    xold = xnew;
    if k > kmax
       my_continue = 1;
       disp('maximum number of iteration reached')
    end
end

function g = grad(x)
    g(1) = 1-(1/x(1)^2);
    g(2) = 1-(1/x(2)^2);
end

function H = hessian(x)
    H = [2/x(1)^3 0;...
         0 2/(x(2)^3)];
end

function f = obj(x)
    f = x(1)+(1/x(1))+x(2)+(1/x(2));
end