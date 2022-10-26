function x_out = scaling(x,l,u,type)
%
% scaling.m scales or unscales the vector x according to the bounds
% specified by u and l. The flag type indicates whether to scale (1) or
% unscale (2) x. Vectors must all have the same dimension.

% ensure same size matrices
ss = size(x);
ss_l = size(l);
if ss(1) == 1 || ss(2) == 1
    if ss(2) ~= ss_l(2)
        msg = 'matrix dimensions incompatible, transposing';
        warning(msg)
        x = x';
    end
end

if type == 1
    % scale
    x_out = (x-l)./(u-l);
elseif type == 2
    % unscale
    x_out = l + x.*(u-l);
end