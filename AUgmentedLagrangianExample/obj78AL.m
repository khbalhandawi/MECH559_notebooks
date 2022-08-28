function f = obj78AL(x,r)
global lamda


f = 3*x(1)*x(1) + x(2)*x(2) + lamda*(2-x(1)-x(2)) + ((2-x(1)-x(2))*(2-x(1)-x(2)))/r; 

return
end