function h = append_point(d,style,color,marker_size,linewidth,fig)
    ss = size(d);
    q = combnk(1:ss(2),2); % choose 2 out d variables
    q = flipud(q);
    ss = size(q);
    for par = 1:1:ss(1)
        % Plot points
        i = q(par,:); % plot variable indices
        ch = fig.Children;
        hold(ch(par*2 + 1),'on');
        h = plot(d(:,i(1)), d(:,i(2)),style,'MarkerEdgeColor',color,'markersize',marker_size,'linewidth',linewidth,'Parent', ch(par*2 + 1));
        uistack(h,'top')
    end
end

