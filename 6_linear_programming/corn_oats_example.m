clear all
clc
close all
format short
format compact

%% Solving the problem by visualization

vertices = [0.,0.; 0.,240.;160.,0.;80.,160.];

% inputs for visualization
n=150;
d=linspace(-40,320,n);
[X1,X2]=meshgrid(d,d);

F=-40*X1-30*X2;

% constraint lines
x = d;
% 2*x1+x2<=320
g1=(320-2*x)/1;
% x1+x2<=240
g2=(240-x)/1;
% x1<=5
g3=x; % vertical line
% x2>=-5
g4=0*ones(size(x)); % horizontal line

% feasible region
G1=2*X1+X2-320;
G2=X1+X2-240;
G3=X1;
G4=X2;

G=(G1<=0).*(G2<=0).*(G3>=0).*(G4>=0);
G=G./G;
F=G.*F;

% initialize figure
f1 = figure(1);
contour(X1,X2,F,21)
hold on;
plot(x,g1,'r',x,g2,'b',0*ones(size(x)),g3,'k',x,g4,'m');
plot(vertices(:,1),vertices(:,2),'.k','MarkerSize',20)
ylim([-40,320])
l = legend('$-40x_1-30x_2$','$2x_1+x_2\leq320$','$x_1+x_2\leq240$','$x_1\geq0$','$x_2\geq0$','vertices','interpreter','latex');
axis equal

%% Using basic solutions
A = [2 1 1 0 ; ...
     1 1 0 1 ];
b = [320 240]'; % remember to make b a column vector
c = [-40 -30 0 0]; % this is c transpose

n = size(A,1);
m = size(A,2);

q = combnk(1:m,n);

ss = size(q);
k = 0;
for i = 1:1:ss(1)
    x_sol = zeros(m,1);
    fprintf('-------------------------------------\n')
    fprintf('Sub matrix for basic variables : [ %s ] \n',sprintf('x%d, ',q(i,:)))
    Ab = A(:,q(i,:)) % Sub matrix of basic solution
    x_sol(q(i,:),1) = Ab\b; % Solve subproblem
    fprintf('x_basic = [ %s ] \n',sprintf('%.2f, ',x_sol(q(i,:),1)))
    fprintf('x = [ %s ] \n',sprintf('%.2f, ',x_sol))
    fprintf('f =  %.2f \n',c*x_sol)

    if any(x_sol<0)
        fprintf('THIS IS AN INFEASIBLE SOLUTION\n')
    end
end

%% Using simplex.m

stop = false;
k = 0;
cf = c;
fprintf("\n==================" + ...
    "\nSimplex iterations\n")
while ~stop
    input("hit ENTER to continue")
    [x_opt,stop,A,b,c] = simplex(A,b,c);

    k = k + 1;
    fprintf('iteration = %d\n',k)
    A, c, b, x_opt
    z_opt = x_opt(1:2)
    f_opt = cf*x_opt
    z_opts(k+1,:) = z_opt;
    fprintf("\n==================\n")
end

plot(z_opts(:,1),z_opts(:,2),'--m','LineWidth',2,'MarkerSize',20,'Marker','*')
figure(f1)