function [f,g] = bb(x)
    P = getappdata(0,'PARAM');
    [f,g,~,~] = Blackbox_call(x,P);
end