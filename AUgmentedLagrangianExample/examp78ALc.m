% MECH 559 - M. Kokkolaras
% McGill University

% Example 7_8 in 2nd edition of Papalambros and Wilde 
% (Augmented Lagrangian)


clear all
clc
clf
global lamda

r = [10 1 .1 .01 .001];
lamda = 0;
x0 = [0 0];

for k = 1:length(r)

x1 = -5:.1:5;
x2 = x1;
for j=1:length(x2)
    g(j) = 2 - x2(j);
    for i = 1:length(x1)
        f(j,i) = 3*x1(i)*x1(i)+x2(j)*x2(j);
        fAL(j,i) = obj78ALc([x1(i),x2(j)]',r(k));
    end
end
%figure(1)
%mesh(x1,x2,f)
%figure(2)
%mesh(x1,x2,fAL)

figure(1)
V = 0:1:20; 
cs = contour(x1,x2,f,V); 
clabel(cs)
hold on
plot(g,x2)

figure(2)
V = 0:1:20; 
cs = contour(x1,x2,fAL,V); 
clabel(cs)
hold on
plot(g,x2)

k
%xstarofr = fminunc(@(x) obj78AL(x,r(k)),x0)
options = optimset('GradObj','on'); 
xstarofr = fminunc(@(x) objgr78ALc(x,r(k)),x0,options)
plot(xstarofr(1),xstarofr(2),'r*')

x0 = xstarofr;
lamda=0;
%lamda = lamda + 2*(2-xstarofr(1)-xstarofr(2))/r(k)

pause
end
