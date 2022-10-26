%*******************************************************************************
% heatshield_Omega:  Omega file for the problem, 'heatshield'.
%*******************************************************************************
function [A,l,u] = bb_Omega(n)

    %Param = getappdata(0,'PARAM');
    
    A = [1 0; 0 1;];
    l=[-5; -5];
    u = [5; 5];
return
