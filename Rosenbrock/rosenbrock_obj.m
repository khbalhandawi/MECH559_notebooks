function f = rosenbrock_obj(x)

f = (1 - x(1))^2 + 100 * (x(2) - x(1)^2)^2;