% MECH 559 - M. Kokkolaras
% McGill University
% Rosenbrock function in two dimensions

clear
%close all
clc

% Plot function

x1=-1.5:.05:1.5;
x2=x1;
for i=1:length(x1)
    for j = 1:length(x2)
        f(i,j) = (1 - x1(j))^2 + 100 * (x2(i) - x1(j)^2)^2;
    end
end

figure(1)
clf
mesh(x1,x2,f)

V = 0:1:20; 
figure(2)
clf
cs = contour(x1,x2,f,V); 
clabel(cs)
hold on
plot(1,1,'r+')

% Choose method (gradient = 1, Newton = 2)
method = 1;

% Ask for intitial guess
%xold = input('Type initial guess as a column vector');
xold  = [-1.2 1]'; % problematic initial guess
%xold = [-1 -1]';
xold = [0 0]';
%xold = [1 0]';
%xold = [1 1.2]';
plot(xold(1),xold(2),'r+')
disp(['objective function value = ',num2str(rosenbrock_obj(xold))])

my_epsilon = .001;
kmax = 5000;

my_continue = 0;
k = 0;
while my_continue == 0
    k = k + 1
    if method == 1
       % Gradient method
       %pick alpha or use exact line search by commenting accordingly
       %alpha = .001;
       alpha = rosenbrock_grad(xold)*rosenbrock_grad(xold)'/ ...
              (rosenbrock_grad(xold)*rosenbrock_hessian(xold)*rosenbrock_grad(xold)');
       xnew = xold - alpha*rosenbrock_grad(xold)'
    else
       % Newton's method
       alpha = 1.;
       xnew = xold - alpha*inv(rosenbrock_hessian(xold))*rosenbrock_grad(xold)'
    end
    disp(['norm of gradient = ', num2str(norm(rosenbrock_grad(xnew)))])
    disp(['objective function value = ',num2str(rosenbrock_obj(xnew))])
    plot(xnew(1),xnew(2),'ro')
    %pause
    if norm(rosenbrock_grad(xnew)) <= my_epsilon
        my_continue = 1;
    end
    xold = xnew;
    if k > kmax
       my_continue = 1;
       disp('maximum number of iteration reached')
   end
end   
plot(xnew(1),xnew(2),'r*')