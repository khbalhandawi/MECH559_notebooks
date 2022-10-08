function [x_opt,stop,A,b,c,k] = simplex(A,b,c,k,verbose)
    %simplex: performs a single simplex iteration
    % INPUTS
    %   A: The coefficient matrix A
    %   b: the solution vector b
    %   c: the objective vector c
    %   k: iteration number
    %   verbose; true/false if true prints the outputs
    % OUTPUTS
    %   x_opt: The optimal variable values
    %   stop: a flag informing the user that no further 
    %       iterations are possible
    %   A: The coefficient matrix A
    %   b: the solution vector b
    %   c: the objective vector c
    %   k: iteration number
    %   verbose; true/false if true prints the outputs

    if nargin == 3
        verbose = true;
    end
    
    ss = size(A);

    % optimal solution
    x_opt = zeros(ss(2),1);
    for j = 1:1:ss(2) % rows
        if c(j) == 0
            x_opt(j) = max(A(:,j)'*b,0);
        end
    end

    if verbose
        fprintf('iteration = %d\n',k)
        A, c, b, x_opt
    end

    [j,q] = min(c); % pivot column
    if j < 0
        A_row_search = b./A(:,q);
        A_row_search(A_row_search<0) = inf;
        [i,p] = min(A_row_search); % pivot row
        if i ~= inf
            An = A; cn = c; bn = b;
            for i = 1:1:ss(1) % rows
                if i == p
                    for j = 1:1:ss(2) % columns
                        A(i,j) = (An(i,j)/An(i,q));
                        b(i) = (bn(i)/An(i,q));
                    end
                else
                    for j = 1:1:ss(2) % columns
                        A(i,j) = An(i,j) - (An(i,q)/An(p,q))*An(p,j);
                        c(j) = cn(j) - (cn(q)/An(p,q))*An(p,j);
                    end
                    b(i) = bn(i) - (An(i,q)/An(p,q))*bn(p);
                end
            end

            stop = false;
            k=k+1;
        else
            fprintf('Problem is unbounded at iteration = %d\n',k)
            x_opt = Inf*ones(ss(2),1);
            stop = true;
        end
    else
        fprintf('solution terminated at iteration = %d\n',k)
        stop = true;
    end
end