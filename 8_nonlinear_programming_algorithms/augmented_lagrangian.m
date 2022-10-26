% MECH 559 - M. Kokkolaras
% McGill University

% Example 7_8 in 2nd edition of Papalambros and Wilde 
% (Augmented Lagrangian method)


clear all
clc
clf
addpath('../../plotting_functions')
addpath('../../plotting_functions/hatchfill2_r8')
global lambda

%% Visualizing the optimization problem
% variables
x1 = linspace(-5,5,100);
x2 = x1;
[X1,X2] = meshgrid(x1,x2);


% objective and constraint
F = 3*X1.*X1+X2.*X2;
G1 = (2 - X2 - X1);

% initialize figure
fig1 = figure(1);
fig1.Position=[100 100 600 500];
ax = axes(fig1);
axis(ax,[-5,5,-5,5]) % fix the axis limits

[cc, h_obj] = contour(ax,X1,X2,F,'LineWidth',2);
hold(ax,'on')
colorbar(ax)
title_obj = ['$f(\mathbf{x},\mathbf{p})$'];

plot_constraint(ax,X1,X2,G1,'r');
cstr_h(1) = plot(ax,NaN,NaN,'linestyle','-','color','r'); % collect contour handles
label_cstr{1} = ['${g_',num2str(1),'}(\mathbf{x};\mathbf{p}) > 0$'];

%% Begin optimization
r = [10 1 .1 .01 .001]; % sequence of decreasing r's
lambda = 0;
x0 = [0 0];
for k = 1:length(r)

    for j=1:length(x2)
        for i = 1:length(x1)
            [fAL(j,i),~] = objgrAL([x1(i),x2(j)]',r(k));
        end
    end
    
    V = 0:1:5; 
    [cc2, h_obj2] = contour(ax,x1,x2,fAL,V,'EdgeColor','k','LineWidth',2); 
    figure(fig1)
    
    k
    options = optimset('GradObj','on'); 
    optfunc = @(x) objgrAL(x,r(k));
    xstarofr = fminunc(optfunc,x0,options)
    h_opt = plot(xstarofr(1),xstarofr(2),'r*','MarkerSize',12,'LineWidth',2);
    
    x0 = xstarofr;
    lambda = lambda + 2*(2-xstarofr(1)-xstarofr(2))/r(k)
    
    pause
end

%% Plot legend
label_opt = ['$\mathbf{x}_k^*$'];
label_obj2 = ['$T(\mathbf{x},r_k)$'];

lh = legend([h_obj cstr_h h_obj2 h_opt],[title_obj label_cstr label_obj2 label_opt],'Orientation','horizontal');
rect = [0.3613, 0.945, .27, .0525];
set(lh, 'Position', rect, 'interpreter', 'latex', 'fontsize', 13)
set(fig1,'color','w');

%% Optimization functions with penalty applied
function [f,g] = objgrAL(x,r)
    global lambda
    f = 3*x(1)^2 + x(2)^2 + lambda*(2-x(1)-x(2)) + max(0,((2-x(1)-x(2))^2)/r); 
    g = [6*x(1)-lambda-max(0,(2-x(1)-x(2))/r)
         2*x(2)-lambda-max(0,(2-x(1)-x(2))/r)];
end