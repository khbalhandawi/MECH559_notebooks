% MECH 559 - M. Kokkolaras
% McGill University
% Conjugate gradients method for Rosenbrock's 2-d function


clear
%close all
clc

% Plot function
x1=-1.5:.05:1.5;
x2=-.5:.05:1.5;
for i=1:length(x2)
    for j = 1:length(x1)
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


% Ask for intitial guess
%xold = input('Type initial guess as a column vector');
xold  = [-1.2 1]'; % problematic initial guess
xold = [-1 -1]';
xold = [0 0]';
%xold = [1 0]';
plot(xold(1),xold(2),'ro')
disp(['objective function value = ',num2str(rosenbrock_obj(xold))])

my_epsilon = .01;
kmax = 200;

my_continue = 0;
k = 1;
d = -rosenbrock_grad(xold)';
while my_continue == 0
      alpha_k = - (rosenbrock_grad(xold)*d)/(d'*rosenbrock_hessian(xold)*d);
      %alpha_k = - (rosenbrock_grad(xold)*d)/(d'*d);
xnew = xold + alpha_k * d;
      k = k + 1
      %beta_k = rosenbrock_grad(xnew) * d / (d'* d);
      beta_k = rosenbrock_grad(xnew) * rosenbrock_hessian(xold) * d / (d'* rosenbrock_hessian(xold) * d);
      %beta_k = rosenbrock_grad(xnew) * rosenbrock_grad(xnew)' / (d' * d);
      %beta_k = rosenbrock_grad(xnew) * d / (d'* rosenbrock_hessian(xnew) * d);
      dnew = - rosenbrock_grad(xnew)' + beta_k * d;
      xold = xnew;
      d = dnew;
      disp(['norm of gradient = ', num2str(norm(rosenbrock_grad(xnew)))])
      disp(['objective function value = ',num2str(rosenbrock_obj(xnew))])
      plot(xnew(1),xnew(2),'ro')
      %pause
      if norm(rosenbrock_grad(xnew)) <= my_epsilon
         disp('gradient norm less than threshold: converged')
         my_continue = 1;
      end
      if k > kmax
         my_continue = 1;
         disp('maximum number of iteration reached')
      end
end  