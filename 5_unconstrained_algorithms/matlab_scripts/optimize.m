function [xnew,fnew] = optimize(obj,grad,hessian,method,xold,my_epsilon, ...
    kmax,interact,verbose,plot_progress)
    % optimize, this function implements:
    %   1) Gradient descent with exact line search
    %   2) Newtons method
    % Inputs:
    %   obj             : a user-defined function for calculating the objective
    %   grad            : a user-defined function for calculating the gradient
    %   hessian         : a user-defined function for calculating the Hessian
    %   method          : (gradient = 1, Newton = 2)
    %   xold            : a vector representing the initial guess
    %   my_epsilon      : stopping tolerance
    %   kmax            : maximum number of iterations allowed
    %   interact        : (on = 1, off = 0), turn this on so that you an interactively see optimization progress by hitting ENTER
    %   verbose         : (on = 1, off = 0), turn this on for printing each iteration
    %   plot_progress   : (on = 1, off = 0), whether to plot optimization progress on current figure
    % Outputs:
    %   xnew        : final point at which termination occurs
    %   fnew        : objective function value at termination
    
    % Default arguments
    if ~exist('my_epsilon','var')
        my_epsilon = 1e-8;
    end
    if ~exist('kmax','var')
        kmax = 1e5;
    end
    if ~exist('interact','var')
        interact = 0;
    end
    if ~exist('verbose','var')
        verbose = 0;
    end
    if ~exist('plot_progress','var')
        plot_progress = 0;
    end

    if plot_progress == 1
        % Get problem dimensionality
        if size(xold,2) == 2 % 2D problems
            plot(xold(1),xold(2),'r*')
        elseif size(xold,2) == 1 % 1D problems
            plot(xold,obj(xold),'r*')
        end
    end

    if verbose == 1
        disp(['objective function value = ',num2str(obj(xold))])
    end
    
    if interact == 1
        pause
    end
    
    my_continue = 0;
    k = 0;
    while my_continue == 0
        k = k + 1;
        if method == 1
           % Gradient method
           alpha = grad(xold)*grad(xold)'/ ...
                   (grad(xold)*hessian(xold)*grad(xold)');
           xnew = xold - alpha*grad(xold)';
        else
           % Newton's method 
           alpha = 1;
           xnew = xold - alpha*inv(hessian(xold))*grad(xold)';
        end
        fnew = obj(xnew);
        
        % plot progress
        if plot_progress == 1
            if size(xold,1) == 2 % 2D problems
                plot(xnew(1),xnew(2),'r*')
            elseif size(xold,1) == 1 % 1D problems
                plot(xnew,obj(xnew),'r*')
            end
        end

        % print progress
        if verbose == 1
            k
            alpha
            xnew
            disp(['norm of gradient = ', num2str(norm(grad(xnew)))])
            disp(['objective function value = ',num2str(fnew)])
        end

        if interact == 1
            pause
        end

        if norm(grad(xnew)) <= my_epsilon
            my_continue = 1;
        end
        xold = xnew;
        if k > kmax
           my_continue = 1;
           disp('maximum number of iteration reached')
       end
    end
end