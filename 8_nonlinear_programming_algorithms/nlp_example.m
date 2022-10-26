addpath('../MATLAB_Algorithms/NomadM')
addpath('../MATLAB_Algorithms/dace')
addpath('../MATLAB_Algorithms')
addpath('../plotting_functions/hatchfill2_r8')
addpath('../plotting_functions')

% Nonlinear programming algorithms
close all
clear all
clc

%% Problem definition
% lower and upper bounds
lb =[-100; -100];
ub =[100; 100];

% initial guess
x0 = [-5; 20];

% matrix of linear inequality constraints
A = [2  3; -5 -2; -2 7];
% right hand side of linear inequality constraints
b = [10; -2; 8];

% matrix of linear equality constraints
Aeq = [];
% right hand side of linear equality constraints
beq = [];

%% Use fmnincon to solve the example problem
my_algorithm = 'sqp';
% my_algorithm = 'interior-point';
% my_algorithm = 'interior-point-exact-gradient';
% my_algorithm = 'direct';

switch my_algorithm
    
    case 'sqp'

    % change algorithm option to SQP 
    my_options = optimoptions('fmincon','Algorithm','sqp','Display','iter',...
        'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);
    % call fmincon using SQP algorithm
    [xopt,fopt,ExitFlag,out,lamda] = fmincon(@my_fun,x0,A,b,Aeq,beq,lb,ub,@my_nonlcon,my_options);
    f_count = out.funcCount;

    case 'interior-point'

    % change algorithm option to interior-point
    my_options = optimoptions('fmincon','Algorithm','interior-point','Display','iter',...
        'SpecifyObjectiveGradient',false,'SpecifyConstraintGradient',false);
    % call fmincon using interior-point algorithm (default)
	[xopt,fopt,ExitFlag,out,lamda] = fmincon(@my_fun,x0,A,b,Aeq,beq,lb,ub,@my_nonlcon,my_options);
    f_count = out.funcCount;

    case 'interior-point-exact-gradient'
        
    %use interior-point with exact derivatives (not finite difference approximations)
    my_options = optimoptions('fmincon','Algorithm','interior-point','Display','iter',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
    [xopt,fopt,ExitFlag,out,lamda] = fmincon(@my_fun,x0,A,b,Aeq,beq,lb,ub,@my_nonlcon,my_options);
    f_count = out.funcCount;

    case 'direct'
        
    % use DIvidedRECTangles
    nc = 1; % number of non-linear constraints 
    c_L = -inf*ones(nc,1); % non-linear constraint lower bounds (= -infinity) 
    c_U = zeros(nc,1); % negative null form
    b_L = -inf*ones(length(b),1); % linear constraint lower bounds (= -infinity)
    % Integer variables
    I=[];
    warning off
    GLOBAL.MaxIter=0;
    GLOBAL.MaxEval=2000;
    %call algorithm
    outputs = gclsolve(@my_fun,@my_nonlcon,lb,ub,A,b_L,b,c_L,c_U,I,GLOBAL); %,PriLev,varargin)
    xopt = outputs.x_k;
    fopt = outputs.f_k;
    f_count = outputs.FuncEv;
    
end

fprintf('================================\n')
fprintf('SOLUTION\n')
fprintf('The optimizer is at x = [%f, %f]\n',xopt(1),xopt(2))
fprintf('The optimum function value is f = %f \n',fopt)
fprintf('The number of function evaluations was = %f \n',f_count)
fprintf('================================\n')

%% plot objective and constraints
fig1 = figure(1);
fig1.Position = [100 100 900 600];

ax = axes(fig1);

resolution = 100;
x1 = linspace(-5,5,resolution);
x2=x1;
[X1,X2] = meshgrid(x1,x2);
F = X1.^2 + (X2-1).^2;

[cf,hf] = contourf(X1,X2,F);
hold(ax,'on')
colorbar(ax)

% plot constraints
G1= -10+2*X1+3*X2;
G2= 2-5*X1-2*X2;
G3= -8-2*X1+7*X2;
G4= -0.2*X1.^2+X2;

plot_constraint(ax,X1,X2,G1,'r');
plot_constraint(ax,X1,X2,G2,'b');
plot_constraint(ax,X1,X2,G3,'g');
plot_constraint(ax,X1,X2,G4,'m');

% create legend labels
ax = gca;
labels = {'$f({x_1},{x_2})$'};
handles = [hf];
colors = {'r','b','g','m'};
for i = 1:1:4
    labels{end+1} = ['$g_',num2str(i),'({x_1},{x_2})>0$'];
    handles(end+1) = plot(ax,NaN,NaN,'linestyle','-','color',colors{i}); % collect contour handles
end

%% plot optimizer
hx = plot(xopt(1),xopt(2),'r*','LineWidth',2,'MarkerSize',12);

labels{end+1} = '$\mathbf{x}^*$';
handles(end+1) = hx;

lh = legend(handles,labels,'Orientation','horizontal');
rect = [0.3613, 0.945, .27, .0525];
set(lh, 'Position', rect, 'interpreter', 'latex', 'fontsize', 13)

%% User defined objective (and gradient) and nonlinear constraints (and their gradients)

function [f,dfdx] = my_fun(x)
    % function to evaluate objective function
    
    f = x(1)^2 + (x(2)-1)^2;
    dfdx = [2*x(1) 2*(x(2)-1)];
    %f=-f;

    % in case of simulation:
    % write x to a text file or as command line arguments
    % call/run simulation (either reads the text file or accepts arguments)
    % read results from a text file (output by the simulation)
    % f = ...
    % see examples in 9_dfo
end

function [gineq,geq,dgineqdx,dgeqdx] = my_nonlcon(x)
    %function that evaluates general nonlinear constraints
    %returns inequality and equality constraints (negative null formulation)
    % and their gradients

    %evaluate constraints and derivatives
    gineq(1) = x(2) - 1/5*x(1)^2;
    geq=[];
    dgineqdx = [-2/5*x(1) 1]';
    dgeqdx = [];
end
