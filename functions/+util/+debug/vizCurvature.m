x_range = min(x):max(x);
y_range = feval(model,out,x_range);
y_range_init = feval(model,pinit,x_range);

figure(figno),
clear plt
if iadj==1
    % range
%     ran = [min(dispvec);max(dispvec)];xl = ran(:,iadj);xl=xl(:)';
%     yl = [0 dims(setdiff([1 2],iadj))];
%     xlim(xl), ylim(yl)
    
    subplot(2,3,[1 3])
    cla
    plt1 = plot(y,x,'+');
    hold on
    plt2 = plot([nan;y(outliers)],[nan;x(outliers)],'ro');
    plt3 = plot(y_range,x_range,'m-','LineWidth',4);
    plt4 = plot(y_range_init,x_range,'g-','LineWidth',2);
    daspect([1 50 1])
    legend([plt1 plt2 plt3 plt4],'matched feats','outliers','estimated model','initial model')
    
elseif iadj==2
    % range
%     ran = [min(dispvec);max(dispvec)];xl = ran(:,iadj);xl=xl(:)';
%     yl = [0 dims(setdiff([1 2],iadj))];
%     xlim(yl), ylim(xl)
    
    subplot(2,3,[2 3])
    cla
    plt1 = plot(x,y,'+');
    hold on
    plt2 = plot([nan;x(outliers)],[nan;y(outliers)],'ro');
    plt3 = plot(x_range,y_range,'m-','LineWidth',4);
    plt4 = plot(x_range,y_range_init,'g-','LineWidth',2);
    daspect([30 1 1])
    legend([plt1 plt2 plt3 plt4],'matched feats','outliers','estimated model','initial model')
end