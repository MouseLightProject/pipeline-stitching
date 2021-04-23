function vector_field_3d_debug_script(scopeloc, scopeparams, params, XYZ_tori, XYZ_tp1ori, outliers)
    Rmin = min(scopeloc.loc)*1e6;
    Rmin=Rmin*.995;
    Rmax = round((max(scopeloc.loc)*1e3+scopeparams(1).imsize_um)*1e3);
    Rmax=Rmax*1.005;
    
    %[bandwidth,density,X,Y]=kde2d(XYZ_tori(:,1:2),256,Rmin(1:2),Rmax(1:2));
    figure(35), clf, cla,
    view([0,90])
    set(gca,'Ydir','reverse')
    hold on
    colormap hot,
    set(gca, 'color', 'w');
    %         plot(XYZ_tori(:,1),XYZ_tori(:,2),'w+','MarkerSize',5)

    plot3(XYZ_tori(:,1),XYZ_tori(:,2),XYZ_tori(:,3),'k.') % layer t
    plot3(XYZ_tp1ori(:,1),XYZ_tp1ori(:,2),XYZ_tp1ori(:,3),'r.') % layer tp1
    myplot3(XYZ_tori(outliers,:),'md') % layer t
    myplot3(XYZ_tori(outliers,:),'gd') % layer t
    x = scopeloc.loc(ix,1)*1e6;
    y = scopeloc.loc(ix,2)*1e6;
    w = scopeparams(1).imsize_um(1)*ones(sum(ix),1)*1e3;
    h = scopeparams(1).imsize_um(2)*ones(sum(ix),1)*1e3;
    for ii=1:sum(ix)
        rectangle('Position', [x(ii)-w(ii) y(ii)-h(ii) w(ii) h(ii)],'EdgeColor','r')
        text(x(ii)-w(ii)/2,y(ii)-h(ii)/2,2,sprintf('%d: %d',ii,idxinlayer(ii)),...
            'Color','m','HorizontalAlignment','center','FontSize',8)
    end

    for ii=find(idxinlayer==5163)
        rectangle('Position', [x(ii)-w(ii) y(ii)-h(ii) w(ii) h(ii)],'EdgeColor','w','LineWidth',2,'FaceColor',[0 .5 .5])
    end
    %surfc(X,Y,density/max(density(:)),'LineStyle','none'),
    set(gca,'Ydir','reverse')
    %         legend('P_t','P_{t+1}','E_t','E_{t+1}')
    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.02 .02]/2 , ...
        'XTickLabel'  , (1:10) , ...
        'YTickLabel'  , (1:10) , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'XGrid'       , 'on'      , ...
        'YGrid'       , 'on'      , ...
        'XColor'      , [.3 .3 .3], ...
        'YColor'      , [.3 .3 .3], ...
        'LineWidth'   , 1         );
    ax = gca;
    ax.YLabel.String = 'mm';
    ax.YLabel.FontSize = 16;
    ax.XLabel.String = 'mm';
    ax.XLabel.FontSize = 16;
    xlim([Rmin(1) Rmax(1)]-params.imsize_um(1)*1e3)
    ylim([Rmin(2) Rmax(2)]-params.imsize_um(2)*1e3)
    % title([num2str(t),' - '])
    drawnow
    %
    vizfolder = 'qualityfold' ;
    if ~exist(vizfolder,'dir')
        mkdir(vizfolder)
    end
    export_fig(fullfile(vizfolder,sprintf('Slice-%05d.png',t)))
end
