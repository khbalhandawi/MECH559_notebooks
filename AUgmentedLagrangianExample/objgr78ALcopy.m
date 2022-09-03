function [f,g] = objgr78AL(x,r)

global lamda

f = 3*x(1)*x(1) + x(2)*x(2) + ((2-x(1)-x(2))*(2-x(1)-x(2)))/r; 
g = [6*x(1)
     2*x(2)];
 
return
end