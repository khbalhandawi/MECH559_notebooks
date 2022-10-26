clearvars
close all
clc
addpath('../MATLAB_Algorithms/NomadM')
addpath('../MATLAB_Algorithms/dace')
addpath('../MATLAB_Algorithms')
addpath('../plotting_functions/hatchfill2_r8')
addpath('../plotting_functions')
addpath('bb_nomadm')

% Keep track of number of function calls and the last called value of x
global index P xLast f_bb gradf_bb g_bb dg_bb

%% Problem definition
parameters = [8,-4,-3,-3,-3,3,2,0.8]; % <---- with these parameters, the constraints are inactive
parameters = [4,-4,-3,-3,-3,3,2,0.8]; % <---- with these parameters, the constraints are active

lb = [-5, -5];
ub = [5, 5];

%% Construct surrogate
train_new_surrogate = false; % run this once only

if train_new_surrogate
    save_data = false;
    use_surrogate = false;
    use_gradients = false;
    sur_model = [];
    P = {save_data,use_surrogate,use_gradients,parameters,sur_model}; 

    n_samples = 50;
    X_train = lhsdesign(n_samples,2);
    X_train = scaling(X_train,lb,ub,2); % unscale from (0,1) to (-5,5)
    for i = 1:1:n_samples
        [f,g,~,~] = Blackbox_call(X_train(i,:),P);
        Y_train(i,:) = [f g'];
    end
    theta = 1*ones(1,2);
    [dmodel, perf] = dacefit(X_train, Y_train, @regpoly0, @correxp , theta);
    save surrogate_model dmodel
else
    load surrogate_model dmodel
end

%% Blackbox optimization problem solved using fmincon
% test = Blackbox_call([1,2],P) % its a good idea to test your blackbox function first
save_data = true;
use_surrogate = false;
use_gradients = false;
sur_model = dmodel;
P = {save_data,false,use_gradients,parameters,sur_model}; 

index = 0; % Initialize counter
xLast = [];
f_bb = [];
gradf_bb = [];
g_bb = [];
dg_bb = [];

x0 = [0, 0];  % <--------------------------------------------------------- SET INITIAL GUESS HERE

% Delete log files (from previous runs)
if save_data
    try_remove('log.txt');
end

% using MATLAB optimization toolbox
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
%     'SpecifyObjectiveGradient',use_gradients,'SpecifyConstraintGradient',use_gradients,...
%     'FiniteDifferenceType','central','PlotFcn',@optimplotfval);
% [x_opt,f_opt,exitflag,output,lambda] = fmincon(@obj,x0,[],[],[],[],lb,ub,@nonlinear_cstr,options);

% % using global optimization toolbox
% options = optimoptions('patternsearch','Display','iter',...
%     'MeshTolerance',1e-3,'ConstraintTolerance',1e-3,...
%     'PlotFcn',@psplotbestf);
% [x_opt,f_opt,exitflag,output] = patternsearch(@obj,x0,[],[],[],[],lb,ub,@nonlinear_cstr,options);

% using NOMADm
[BestF,BestI,RunStats,RunSet] = mads_batch;
x_opt = BestF.x;
f_opt = BestF.f;

fprintf('================================\n')
fprintf('SOLUTION\n')
fprintf('The optimizer is at x = [%f, %f]\n',x_opt(1),x_opt(2))
fprintf('The optimum function value is f = %f \n',f_opt)
fprintf('The number of blackbox evaluations was = %f \n',index)
fprintf('================================\n')

%% Visualize results
n_divisions = 50;

x = linspace(lb(1),ub(1),n_divisions);
y = linspace(lb(2),ub(2),n_divisions);

[X1,X2] = meshgrid(x,y);

% evaluate constraints and objective
if use_surrogate
    X = [reshape(X1,[],1) reshape(X2,[],1)];
    [YX, MSE] = predictor(X, dmodel);

    F = reshape(YX(:,1),size(X1));
    G1 = reshape(YX(:,2),size(X1));
    G2 = reshape(YX(:,3),size(X1));
    G3 = reshape(YX(:,4),size(X1));
else
    obj_fun = @(x,y) (((x.^2) + y - 11) .^2) + ((x + (y.^2) -7) .^2);
    F = obj_fun(X1,X2);
    G1 = constraint_fun1([reshape(X1,[],1),reshape(X2,[],1)]);
    G2 = constraint_fun2([reshape(X1,[],1),reshape(X2,[],1)]);
    G3 = constraint_fun3([reshape(X1,[],1),reshape(X2,[],1)]);
    G1 = reshape(G1,size(X1));
    G2 = reshape(G2,size(X1));
    G3 = reshape(G3,size(X1));
end

% plot the function
fig1 = figure(2);
fig1.Position=[100 100 900 600];
ax = axes(fig1);
axis(ax,[lb(1),ub(1),lb(2),ub(2)]) % fix the axis limits

[cc, h_obj] = contourf(ax,X1,X2,F);
hold(ax,'on')
colorbar(ax)
title_obj = ['$f(\mathbf{x},\mathbf{p})$'];

plot_constraint(ax,X1,X2,G1,'r');
cstr_h(1) = plot(ax,NaN,NaN,'linestyle','-','color','r'); % collect contour handles
label_cstr{1} = ['$\hat{g_',num2str(1),'}(\mathbf{x};\mathbf{p}) > 0$'];

plot_constraint(ax,X1,X2,G2,'c');
cstr_h(2) = plot(ax,NaN,NaN,'linestyle','-','color','c'); % collect contour handles
label_cstr{2} = ['$\hat{g_',num2str(2),'}(\mathbf{x};\mathbf{p}) > 0$'];

plot_constraint(ax,X1,X2,G3,'g');
cstr_h(3) = plot(ax,NaN,NaN,'linestyle','-','color','g'); % collect contour handles
label_cstr{3} = ['$\hat{g_',num2str(3),'}(\mathbf{x};\mathbf{p}) > 0$'];

lh = legend([h_obj cstr_h],[title_obj label_cstr],'Orientation','horizontal');
rect = [0.3613, 0.945, .27, .0525];
set(lh, 'Position', rect, 'interpreter', 'latex', 'fontsize', 13)

set(fig1,'color','w');
% export_fig('bb_opt.pdf','-p0.002',fig1); 

%% Visualize optimization
cyan = [0 1 1]; yellow = [1 1 0];

plot(ax,x0(1),x0(2),'x','MarkerEdgeColor','c','markersize',8,'linewidth',1.5);
plot(ax,[x0(1),x_opt(1)],[x0(2),x_opt(2)],'-x','MarkerEdgeColor','y','markersize',8,'linewidth',1.5);

lh = legend([h_obj cstr_h],[title_obj label_cstr],'Orientation','horizontal');
rect = [0.3613, 0.945, .27, .0525];
set(lh, 'Position', rect, 'interpreter', 'latex', 'fontsize', 13)

set(fig1,'color','w');
% export_fig('5h_opt_result_2.pdf','-p0.002',fig1); 

%% Objective and constraint functions

% These functions are much slower to use because the solver calls the black
% box twice! once for the objective and once for the constraints
% function [f,gradf] = obj(x)
%     global P
%     [f,~,gradf,~] = Blackbox_call(x,P);
% end
% 
% function [cineq,ceq,Dcineq,Dceq] = nonlinear_cstr(x)
%     global P
%     [~,g,~,dg] = Blackbox_call(x,P);
%     cineq = g;
%     Dcineq = dg';
%     ceq = [];
%     Dceq = [];
% end

% These functions are faster to call
function [f,gradf] = obj(x)
    global P xLast f_bb g_bb gradf_bb dg_bb
    if ~isequal(x,xLast) % Check if computation is necessary
        [f_bb,g_bb,gradf_bb,dg_bb] = Blackbox_call(x,P);
        xLast = x;
    end
    % Now compute objective function
    f = f_bb;
    gradf = gradf_bb;
end

function [cineq,ceq,Dcineq,Dceq] = nonlinear_cstr(x)
    global P xLast f_bb g_bb gradf_bb  dg_bb
    if ~isequal(x,xLast) % Check if computation is necessary
        [f_bb,g_bb,gradf_bb,dg_bb] = Blackbox_call(x,P);
        xLast = x;
    end
    % Now compute constraint function
    cineq = g_bb;
    Dcineq = dg_bb';
    ceq = [];
    Dceq = [];
end

% These functions are for plotting only
function [cstr] = constraint_fun1(x)
    global P
    p = P{4};
    h = p(1); d = p(2);
    cstr = - ( (x(:,1) - h).^2 - (x(:,2) - d) - 5);    
end

function [cstr] = constraint_fun2(x)
    global P
    p = P{4};
    h = p(3); d = p(4);
    cstr = -( ((x(:,1) - h).^2) + ((x(:,2) - d).^2) - 4 );
end

function [cstr] = constraint_fun3(x)
    global P
    p = P{4};
    h = p(5); d = p(6);
    r1 = p(7); r2 = p(8);
    cstr = -( ((r1*(x(:,1) - h)).^2) + ((r2*(x(:,2) - d)).^2) - 8 );
end

function try_remove(filename)
    % Delete files
    if exist(filename, "file") == 2
        delete(filename)
    end
end