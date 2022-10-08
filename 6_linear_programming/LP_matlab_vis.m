% MECH 559 
% McGill University; M. Kokkolaras
% Linear programming example

clear all
close all
clc

% Example problem
A = [...
    2 3 -2 -3 1 0 0 0 0; ...
    -5 -2 5 2 0 1 0 0 0; ...
    -2 7 2 -7 0 0 1 0 0; ...
    0 1 0 -1 0 0 0 -1 0; ...
    1 0 -1 0 0 0 0 0 1];
b = [10 -2 8 -5 5]'; % remember to make b a column vector
c = [-1 1 1 -1 0 0 0 0 0]; % this is c transpose

% inputs for visualization
n=150;
x=linspace(-1,6,n);
y=linspace(-6,6,n);

[X,Y]=meshgrid(x,y);

F=c(1)*Y+c(2)*X;

% constraint lines
g1= (10-2*x)./3;
g2= (2-5*x)./2;
g3= (8+2*x)./7;
g4=y;
g5=-5*ones(size(x));

% feasible region
G1=2*X+3*Y-10;
G2=-5*X-2*Y+2;
G3=-2*X+7*Y-8;
G4=X-5;
G5=-Y-5;

G=(G1<=0).*(G2<=0).*(G3<=0).*(G4<=0).*(G5<=0);
G=G./G;
F=G.*F;

% initialize figure
f1 = figure(1);
contour(X,Y,F,21)
hold on;
plot(x,g1,'r',x,g2,'b',x,g3,'g',5*ones(size(x)),g4,'k',x,g5,'m');
ylim([-6,6])
axis equal
input("hit ENTER to continue")
figure(f1) % bring figure window to top

% simplex solution in standard form 
stop = false;
k = 0;
cf = c;
fprintf("\n==================" + ...
    "\nSimplex iterations\n")
while ~stop
    input("hit ENTER to continue")
    [x_opt,stop,A,b,c,k] = simplex(A,b,c,k,true);
    z_opt = x_opt(1:2) - x_opt(3:4)
    f_opt = cf*x_opt
    plot(z_opt(1),z_opt(2),'*r');
    fprintf("\n==================\n")
    figure(f1) % bring figure window to top
end

