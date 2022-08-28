function g = rosenbrock_grad(x)

g(1) = -2*(1 - x(1)) - 400*(x(2) - x(1)^2)*x(1);
g(2) = 200*(x(2) - x(1)^2);