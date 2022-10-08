clear all
clc
close all
format short
format compact

% A = [2 3 1 0 0;-5 2 0 1 0; -2 7 0 0 1];
% b = [ 10 -2 8]';
% c = [ -1 1 0 0 0 ];

% A = [2 1 1 1 0 0;4 2 3 0 1 0 ; 2 5 5 0 0 1];
% b = [ 14 28 30 ]';
% c = [ -1 -2 1 0 0 0];

% A = [2 1 1 1 0 0;1 2 3 0 1 0 ; 2 2 1 0 0 1];
% b = [ 2 5 6 ]';
% c = [ 3 -1 -3 0 0 0];

% % corn/oats example
% A = [1 1 1 0; 2 1 0 1];
% b = [ 240 320 ]';
% c = [ -40 -30 0 0];

% kauserwise example
A = [10 20 1 0; 8 8 0 1];
b = [ 120 80 ]';
c = [ -12 -16 0 0];

% Kokkolaras example
A = [2  3 1 0 0; -5 -2 0 1 0; -2 7 0 0 1];
b = [10; -2; 8];
c = [1 -1 0 0 0];

% A = [-5 -2 1 0; -2 7 0 1];
% b = [-2; 8];
% c = [1 -1 0 0];

xopt = linprog(-c,[],[],A,b)

% Alin = A(:,1:1:3)
% Clin = c(1:1:3)

ss = size(A);
fprintf('iteration = %d\n',0)
A, c, b

for k = 1:1:1000
    [j,q] = min(c); % pivot column
    if j < 0
        A_row_search = b./A(:,q);
        A_row_search(A_row_search<0) = inf
        [i,p] = min(A_row_search); % pivot row
        if i ~= inf
            A_n = [A];
            c_n = c;
            b_n = b;
            for i = 1:1:ss(1) % columns
                if i == p
                    for j = 1:1:ss(2) % rows
                        A_n(i,j) = (A(i,j)/A(i,q));
                        b_n(i) = (b(i)/A(i,q));
                    end
                else
                    for j = 1:1:ss(2) % rows
                        A_n(i,j) = A(i,j) - (A(i,q)/A(p,q))*A(p,j);
                        c_n(j) = c(j) - (c(q)/A(p,q))*A(p,j);
                    end
                    b_n(i) = b(i) - (A(i,q)/A(p,q))*b(p);
                end
            end
            A = A_n; c = c_n; b = b_n;
            fprintf('iteration = %d\n',k)
            A, c, b
            pause
        else
            fprintf('Problem is unbounded at iteration = %d\n',k-1)
            break
        end
    else
        fprintf('solution terminated at iteration = %d\n',k-1)
        break
    end
end