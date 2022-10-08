% MECH 559 
% McGill University; M. Kokkolaras
% Linear programming example

clear all
close all
clc

% inputs for visualization
c = [1 -1];
A = [2  3; -5 -2; -2 7];
b = [10; -2; 8];

n=150;
n1=5;
x=linspace(-n1,n1,n);
y=linspace(-n1,n1,n);

size(x);
[X,Y]=meshgrid(x,y);

F=c(1)*Y+c(2)*X;

g1= (10-2*x)./3;
g2= (2-5*x)./2;
g3= (8+2*x)./7;
g3a = (2*x+8)./7;
g4=0;

G1=2*X+3*Y-10;
G2=-5*X-2*Y+2;
G3=-2*X+7*Y-8;

G=(G1<=0).*(G2<=0).*(G3<=0);
G=G./G;
F=G.*F;

% initialize figure
figure(1)
contour(X,Y,F,21)
hold on;
plot(x,g1,'r', x,g2,'b', x, g3,'g');
pause

% simplex solution in standard form 
% (negative values allowed)
A = [2 3 -2 -3 1 0 0; -5 -2 5 2 0 1 0; -2 7 2 -7 0 0 1];
b = [10; -2; 8];
c = [1 -1 -1 1 0 0 0];

% with x2 >= -5 and x1 <= 5
A = [...
    2 3 -2 -3 1 0 0 0 0; ...
    -5 -2 5 2 0 1 0 0 0; ...
    -2 7 2 -7 0 0 1 0 0; ...
    0 1 0 -1 0 0 0 -1 0; ...
    1 0 -1 0 0 0 0 0 1];
b = [10; -2; 8; -5; 5];
c = [-1 1 1 -1 0 0 0 0 0];

stop = false;
k = 0;
cf = c;
fprintf("\n==================" + ...
    "\nSimplex iterations\n")
while ~stop
    [x_opt,stop,A,b,c,k] = simplex(A,b,c,k,true);
    f_opt = cf*x_opt
    z_opt = x_opt(1:2) - x_opt(3:4);
    plot(z_opt(1),z_opt(2),'*r');
    fprintf("\n==================\n")
    pause
end

