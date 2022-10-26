% MECH 559 - M. Kokkolaras
% McGill University
% Plot reduced functions of examples 5.8 and 5.10 in POD

clear all
clc
close all
clf

x2 = -2:.1:2;
x3 = -3:.1:3;

[X,Y] = meshgrid(x2,x3);
f = -5*Y - 4*X.^2;
F = -5*Y - 4*X.^2;
F = X.^2 + 2*X.*Y;
f=F;

x2star = -95/588;
x3star = 5/14;
fstar = -5*x2star-4*x3star^2;

g = -(x3+x3.^2)./3;
g = 3/2 - x3/2;
G1 = 3*Y + X + X.^2;
G1 = 2*X + Y -3;
G = (G1<=0);
Geq = (G1==0);
%G=G./G;
F = G.*F;
 
figure(1)
mesh(x2,x3,f)
hold on
mesh(x2,x3,G1);

figure(2)
clf
V = -25:1:25; 
cs = contour(x2,x3,f,V); 
clabel(cs)
hold on
plot(x3,g,'r')
plot(x3star,x2star,'*r')

ff = -(3*X + Y - X.^2 - X.*Y - Y.*Y);
figure(3)
mesh(x2,x3,ff)


figure(4)
contour(X,Y,F,21)
hold on
plot(x3,g)
plot(x3star,x2star,'*r')

figure(5)
ff = 5*x3-7*x3.^2;
ffstar = 5*x3star-7*x3star^2;
plot(x3,ff)
hold on
plot(x3star,ffstar,'*r')



Bordered_Hessian = [0 1 1 1; 1 0 -1 -1; 1 -1 0 -1; 1 -1 -1 0]
det(Bordered_Hessian)
stlpm=[0 1 1; 1 0 -1; 1 -1 0];
det(stlpm)
