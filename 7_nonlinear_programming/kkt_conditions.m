clc
clearvars
close all
format compact

%% Class exercise 4 (Q1)
fprintf("====================================\n")
fprintf("Class exercise 4 (Q1):\n\n")
syms x1 x2 m1 m2

f = -x1;
g1 = x2 - (1-x1)^3;
g2 = -x2;
L = f + m1*g1 + m2*g2;

dL = [diff(L,x1);diff(L,x2)];

eq1 = diff(L,x1);
eq2 = diff(L,x2);
eq3 = m1*g1;
eq4 = m2*g2;

S=solve(eq1,eq2,eq3,eq4,x1,x2,m1,m2);

fprintf("the number of solutions is %i\n",length(S.x1))

% we can plot the problem to investigate

x1 = linspace(0,1.7,40);
x2 = linspace(-1.1,1.1,40);

[X1,X2] = meshgrid(x1,x2);

f = -X1;

G1 = (1-x1).^3;
G2 = zeros(1,length(x1));

[C,h] = contourf(X1,X2,f,10);
hold on
p1 = plot(x1,G1,'-r','linewidth',2);
p2 = plot(x1,G2,'-b','linewidth',2);
p3 = plot(1,0,'*m','markersize',10,'linewidth',2);

ax = gca;
labels = {'$f({x_1},{x_2})$',...
    '$g_1({x_1},{x_2})=0$',...
    '$g_2({x_1},{x_2})=0$',...
    '$\mathbf{x}^*$'};
handles = [h;p1;p2;p3];

l = legend(ax,handles, labels,'location','northeast');
set(l,'Interpreter','Latex','fontsize',14);

xlabel('$x_1$','interpreter','LaTex','fontsize',14)
ylabel('$x_2$','interpreter','LaTex','fontsize',14)
xlim([x1(1),x1(end)])
ylim([x2(1),x2(end)])
colorbar

%% Class exercise 3 (Q1)
fprintf("====================================\n")
fprintf("Class exercise 3 (Q1):\n\n")

syms x1 x2 x3 m l

f = -3*x1 + x2 - x3^2;
g = x1 + x2 + x3;
h = -x1 + 2*x2 + x3^2;
L = f + m*g + l*h;

dL = [diff(L,x1);diff(L,x2);diff(L,x3)];

eq1 = diff(L,x1);
eq2 = diff(L,x2);
eq3 = diff(L,x3);
eq4 = l*h;
eq5 = m*g;

S=solve(eq1,eq2,eq3,eq4,eq5,x1,x2,x3,l,m);

fprintf("the number of solutions is %i\n",length(S.x1))

for i=1:1:length(S.x1)
    fprintf("Solution %i: x1 = %.2f, x2 = %.2f, x3 = %.2f, l = %.2f, m = %.2f",...
        i,S.x1(i),S.x2(i),S.x3(i),S.l(i),S.m(i))
    if S.m(i) < 0
        fprintf(", INVALID MULTIPLIERS\n")
    elseif subs(g,[x1,x2,x3],[S.x1(i),S.x2(i),S.x3(i)]) > 0
        fprintf(", VIOLATES CONSTRAINTS\n")
    else
        fprintf("\n")
    end
end

% calculate the Hessian
dL2 = [diff(dL(1),x1),diff(dL(1),x2),diff(dL(1),x3);...
       diff(dL(2),x1),diff(dL(2),x2),diff(dL(2),x3);...
       diff(dL(3),x1),diff(dL(3),x2),diff(dL(3),x3)]

%% Class exercise 4 (Q2)
fprintf("====================================\n")
fprintf("Class exercise 4 (Q2):\nUsing the reduced gradient method:\n\n")

% using the reduced gradients method
syms x1 x2 x3

f = x1^2 + x2^2 + x3^2;
h1 = (x1^2)/4 + (x2^2)/5 + (x3^2)/25 - 1;
h2 = x1 + x2 - x3;

df = [diff(f,x1);diff(f,x2);diff(f,x3)];


h = [h1;h2];
dh = [diff(h(1),x1),diff(h(1),x2),diff(h(1),x3);...
       diff(h(2),x1),diff(h(2),x2),diff(h(2),x3)];

% choose state and decision variables
dhs = dh(:,1:2);
dhd = dh(:,3:end);
dfs = df(1:2,:);
dfd = df(3:end,:);

dw = dfd.' - dfs.'*dhs^(-1)*dhd

eq1 = dw;
eq2 = h1;
eq3 = h2;

S=solve(eq1,eq2,eq3,x1,x2,x3);

fprintf("the number of solutions is %i\n",length(S.x1))

for i=1:1:length(S.x1)
    fprintf("Solution %i: x1 = %.2f, x2 = %.2f, x3 = %.2f\n",...
        i,S.x1(i),S.x2(i),S.x3(i))
    fprintf("Objective: f = %.2f\n",...
        subs(f,[x1,x2,x3],[S.x1(i),S.x2(i),S.x3(i)]))
end

% (This is extra just for your understanding) solving the problem using KKT conditions
fprintf("------------------------------------\n")
fprintf("Using the KKT conditions:\n\n")

syms x1 x2 x3 l1 l2

f = x1^2 + x2^2 + x3^2;
h1 = (x1^2)/4 + (x2^2)/5 + (x3^2)/25 - 1;
h2 = x1 + x2 - x3;
L = f + l1*h1 + l2*h2;

dL = [diff(L,x1);diff(L,x2);diff(L,x3)];
% calculate the Hessian
dL2 = [diff(dL(1),x1),diff(dL(1),x2),diff(dL(1),x3);...
       diff(dL(2),x1),diff(dL(2),x2),diff(dL(2),x3);...
       diff(dL(3),x1),diff(dL(3),x2),diff(dL(3),x3)];

eq1 = diff(L,x1);
eq2 = diff(L,x2);
eq3 = diff(L,x3);
eq4 = l1*h1;
eq5 = l2*h2;

S=solve(eq1,eq2,eq3,eq4,eq5,x1,x2,x3,l1,l2);

fprintf("the number of solutions is %i\n",length(S.x1))

for i=1:1:length(S.x1)
    fprintf("Solution %i: x1 = %.2f, x2 = %.2f, x3 = %.2f, l1 = %.2f, l2 = %.2f",...
        i,S.x1(i),S.x2(i),S.x3(i),S.l1(i),S.l2(i))
    if any([S.l1(i) == 0, S.l2(i) == 0])
        fprintf(", INVALID MULTIPLIERS\n")
    else
        fprintf("\nObjective: f = %.2f\n",...
            subs(f,[x1,x2,x3],[S.x1(i),S.x2(i),S.x3(i)]))
        fprintf("Hession:\n")
        dL2_x_opt = subs(dL2,[x1,x2,x3,l1,l2],[S.x1(i),S.x2(i),S.x3(i),S.l1(i),S.l2(i)])
    end
end