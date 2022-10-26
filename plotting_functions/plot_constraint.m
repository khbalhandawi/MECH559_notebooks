function YX_cstr = plot_constraint(ax,X,Y,YX_cstr,color)
    % Function to plot constraint as a hatched region
    cond = (YX_cstr > 0); % replace infeasible region with zeros
    YX_cstr(cond) = 0;
    if all(YX_cstr == 0,'all')
        YX_cstr(1,1) = -0.1; % override controur not rendered warning
    end

    [c1, h1] = contourf(ax,X,Y,YX_cstr,[0 0]);
    % if not(all(YX_cstr == 0,'all'))
    if not(all(YX_cstr < 0,'all')) % avoid not rendering when entire space is feasible
        set([h1],'linestyle','-','LineColor',color,'LineWidth',2,'Visible','on','Tag','HatchingRegion');
    end

    try
        hatchfill2(h1,'single','LineWidth',2,'Fill','off','HatchVisible','on');
    catch
        fprintf('No area to hatch in constraint\n')
    end
end

