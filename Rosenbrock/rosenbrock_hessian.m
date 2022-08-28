function H = rosenbrock_hessian(x)

H(1,1) = 2 - 400*(x(2) - x(1)^2) + 800*x(1)^2;
H(1,2) = -400*x(1);
H(2,1) = -400*x(1);
H(2,2) = 200;