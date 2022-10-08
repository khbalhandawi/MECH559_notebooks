% MECH 559 
% McGill University; M. Kokkolaras
% Linear programming example

clear all
close all

% vector of objective function coefficients
% c = [-1 1]; % not bounded (bottom right)
c = [1 -1]; % bounded (top left)
%c = [-4 8/7]; % parallel
%Matrix of linear inequality constraints
A = [2  3; -5 -2; -2 7];
% A = [2  3; -5 -2; 2 -7];
%right hand side of linear inequality constraints
b = [10; -2; 8];
%b = [10; -2; -8];

% my_options = optimoptions('linprog');
my_options = optimoptions('linprog','Display','iter');
%my_options = optimoptions('linprog','Algorithm','interior-point','Display','iter');

xopt = linprog(c,A,b)

%pause

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
%G3a = 2*X-7*Y+8;

G=(G1<=0).*(G2<=0).*(G3<=0);
%G=(G1<=0).*(G2<=0).*(G3a<=0);
G=G./G;
F=G.*F;


contour(X,Y,F,21)

hold on;

plot(x,g1,'r', x,g2,'b', x, g3,'g');

