function commit_hatch(h)
    % apply hatch to any tagged regions in plot_constraint
    hp = findobj(h,'Tag','HatchingRegion');
    try
        hatchfill2(hp,'single','LineWidth',2,'Fill','off','HatchVisible','on');
    catch
        fprintf('No area to hatch in constraint\n')
    end
end

