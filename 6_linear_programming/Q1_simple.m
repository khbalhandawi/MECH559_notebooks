clearvars
clc
format compact
close all

A = [2 3 1 0 0;-5 -2 0 1 0; -2 7 0 0 1];
b = [ 10 -2 8]';
c = [ -1 1 0 0 0 ];

q = combnk(1:5,3);

ss = size(q);
k = 0;
for n = 1:1:ss(1)
    xn = zeros(5,1);
    fprintf('-------------------------------------\n')
    A
    fprintf('Sub matrix for basic variables: x%d, x%d, x%d \n',q(n,1),q(n,2),q(n,3))
    Ax = A(:,q(n,:)) % Sub matrix of basic solution
    det(Ax)
    xn(q(n,:),1) = Ax\b; % Solve subproblem
    fprintf('[ x%d, x%d, x%d ] = [ %f, %f, %f ] \n',q(n,1),q(n,2),q(n,3),...
        xn(q(n,1),1),xn(q(n,2),1),xn(q(n,3),1))
    fn = -xn(1) + xn(2) + 0*(xn(3) + xn(4) + xn(5)); % Objective
    % Contraints
    g1(n) = 2*xn(1) + 3*xn(2) -10;
    g2(n) = -5*xn(1) - 2*xn(2) + 2;
    g3(n) = -2*xn(1) + 7*xn(2) - 8;
    
    if g1(n) <= 0 & g2(n) <= 0 & g3(n) <= 0
        k = k+1;
         
        x(:,k) = xn; % Solve subproblem
        f(k) = fn; % Objective
        fprintf('function value %f\n',f(k))
        fprintf('g1 = %f, g2 = %f, g3 = %f\n',g1(n),g2(n),g3(n))
    else
        fprintf('Infeasible solution\n')
        fprintf('g1 = %f, g2 = %f, g3 = %f\n',g1(n),g2(n),g3(n))
    end
end
fprintf('-------------------------------------\n')
([g1',g2',g3']); % Print constraints

[minf,i] = min(f);
fprintf('minimimum function value: f(x*) = %f\n',f(i))
fprintf('minimizer: x* = [%f , %f, %f, %f , %f]\n',x(1,i),x(2,i),x(3,i),x(4,i),x(5,i))

