addpath('../Algorithms/NomadM')
addpath('../Algorithms')

% Nonlinear programming algorithms

close all
clear all
clc

% plot objective and constraints
x1 = -5:.1:5;
x2=x1;
for i=1:length(x1)
    for j = 1:length(x2)
        f(i,j) = x1(i)^2+(x2(j)-1)^2;
    end
end
V = 0:1:20; 
cs = contour(x1,x2,f,V); 
clabel(cs)
hold on
g1= (10-2*x1)./3;
g2= (2-5*x1)./2;
g3= (8+2*x1)./7;
g4= 0.2*x1.*x1;
plot(x1,g1,'r', x1,g2,'b', x1, g3,'g',x1,g4,'m');

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

% lower and upper bounds
lb =[-100; -100];
ub =[100; 100];

% Use fmnincon to solve the example problem

% get current fmincon options 
current_options = optimoptions('fmincon');
% display every iteration
my_options = optimoptions('fmincon','Display','iter');


my_algorithm = 'sqp';
my_algorithm = 'interior-point';
my_algorithm = 'IPED';
my_algorithm = 'direct';

switch my_algorithm
    
    case 'sqp'

    % change algorithm option to SQP 
    my_options = optimoptions('fmincon','Algorithm','sqp','Display','iter');
    % call fmincon using SQP algorithm
    [xoptSQP,foptSQP,ExitFlagSQP,outSQP,lamdaSQP] = fmincon('my_fun',x0,A,b,Aeq,beq,lb,ub,'my_nonlcon',my_options);
    % plot minimizer
    plot(xoptSQP(1),xoptSQP(2),'r*')

    case 'interior-point'

    % change algorithm option to interior-point
    my_options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
    % call fmincon using interior-point algorithm (default)
	[xoptIP,foptIP,ExitFlagIP,outIP,lamdaIP] = fmincon('my_fun',x0,A,b,Aeq,beq,lb,ub,'my_nonlcon',my_options);
    % plot minimizer
    plot(xoptIP(1),xoptIP(2),'r*')
    
    case 'IPED'
        
    %use interior-point with exact derivatives (not finite difference approximations)
    my_options = optimoptions('fmincon','Algorithm','interior-point','Display','iter','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
    [xoptIPED,foptIPED,ExitFlagIPED,outIPED,lamdaIPED] = fmincon('my_funED',x0,A,b,Aeq,beq,lb,ub,'my_nonlconED',my_options);
    % plot minimizer
    plot(xoptIPED(1),xoptIPED(2),'r*')
    
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
    xopt_direct = gclsolve('my_fun','my_nonlcon',lb,ub,A,b_L,b,c_L,c_U,I,GLOBAL); %,PriLev,varargin)
    % plot minimizer
    plot(xopt_direct.x_k(1),xopt_direct.x_k(2),'r*')
    
end


