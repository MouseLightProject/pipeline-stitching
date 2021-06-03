function vizCurvature(x, y, model, out, pinit, outliers, iadj)
    persistent fig ax
    
    if isempty(fig) || ~isvalid(fig) ,
        fig = figure('Color', 'w', 'Name', mfilename()) ;
    end
    if isempty(ax) || ~isvalid(ax) ,
        ax = axes(fig) ;
    end
    cla(ax) ;
    
    x_range = min(x):max(x);
    y_range = feval(model,out,x_range);
    y_range_init = feval(model,pinit,x_range);

    if iadj==1
        plt1 = plot(ax, y,x,'+');
        hold(ax,'on') ;
        plt2 = plot(ax, [nan;y(outliers)],[nan;x(outliers)],'ro');
        plt3 = plot(ax, y_range,x_range,'m-','LineWidth',4);
        plt4 = plot(ax, y_range_init,x_range,'g-','LineWidth',2);
        %daspect([1 50 1])
        hold(ax,'off') ;
        legend(ax, [plt1 plt2 plt3 plt4],'matched feats','outliers','estimated model','initial model')
    elseif iadj==2
        plt1 = plot(ax, x,y,'+');
        hold(ax,'on') ;
        plt2 = plot(ax, [nan;x(outliers)],[nan;y(outliers)],'ro');
        plt3 = plot(ax, x_range,y_range,'m-','LineWidth',4);
        plt4 = plot(ax, x_range,y_range_init,'g-','LineWidth',2);
        hold(ax,'off') ;
        %daspect([30 1 1])
        legend(ax, [plt1 plt2 plt3 plt4],'matched feats','outliers','estimated model','initial model')
    end
    drawnow
end
