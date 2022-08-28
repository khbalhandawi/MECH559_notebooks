function g = ex1grad(x)

g(1) = 4*x(1)^3 - 4*x(1)*x(2);
g(2) = 2*x(2) - 2*x(1)^2;