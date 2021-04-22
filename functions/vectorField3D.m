function vecfield = vectorField3D(params,scopeloc,regpts,scopeparams,curvemodel,theselayers)
    fprintf('Calculating vector fields\n') ;
    Nlayer = params.Nlayer;
    Npts = (params.Ndivs+1).^2;
    dims = params.imagesize;
    order = params.order;
    zlimdefaults = [2 20 dims(3)-21 dims(3)-3];
    params.zlimdefaults = zlimdefaults;
    %beadparams = params.beadparams;
    params.applyFC=1;
    
    tileneighbors = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
    % checkthese = [1 4 5 7]; % 0 - right - bottom - below
    % tileneighbors = tileneighbors(:,checkthese);
    
    % initialize points
    tile_count = length(scopeloc.filepath) ;
    % xy crop for minimal overlap
    [st,ed] = util.getcontolpixlocations(scopeloc,params,scopeparams) ;
    [corrctrlpnttmp_,xlim_cntrl,ylim_cntrl] = util.getCtrlPts(params.imagesize(1),params.imagesize(2),params,st,ed) ;
    % zlimdefaults = [5 25 dims(3)-26 dims(3)-6];
    subcorrctrlpnttmp = corrctrlpnttmp_ + 1 ;
    
    afftransform = reshape([scopeparams(:).affineglFC],3,3,[]);
    % flip x/y dimensions to match imaging
    mflip = -eye(3);%mflip(3,3)=1;
    for ii = 1 : size(afftransform,3) ,
        afftransform(:,:,ii) = afftransform(:,:,ii) * mflip ;
    end

    % for efficiency reasons precompute this
    xlocs = 1:dims(1);
    ylocs = 1:dims(2);
    [xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
    xy = [xy1(:),xy2(:)];

    %
    control = zeros(Npts*Nlayer, 3, tile_count);
    topcntrl1 = zeros(Npts, 3, tile_count);
    topcntrl2 = zeros(Npts, 3, tile_count);
    botcntrl1 = zeros(Npts, 3, tile_count);
    botcntrl2 = zeros(Npts, 3, tile_count);
    zlim_cntrl = zeros(Nlayer, tile_count);
    afftile = zeros(3, 4, tile_count);
    pixstats = nan(tile_count,8);

    filep = strsplit(scopeloc.filepath{1},'/');
    vecfield.root = fullfile('/',filep{1:end-4});

    % INIT
    for tile_index = 1 : tile_count ,
        % copy path
        filepath = fileparts(scopeloc.filepath{tile_index});
        vecfield.path{tile_index} = filepath(length(vecfield.root)+1:end);
        % form an affine matrix
        tform = [afftransform(:,:,tile_index) scopeloc.loc(tile_index,:)'*1e6];
        afftile(:,:,tile_index) = tform;

        % field curvature
        %[locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,tile_index),order,xy,dims,subcorrctrlpnttmp);
        locs = util.fcshift(curvemodel(:,:,tile_index),order,xy,dims,subcorrctrlpnttmp);
        corrctrlpnttmp = locs-1;

        % top layer at z=15: as we have some noisy data at the top
        zlim_1 = zlimdefaults(1);
        corrctrlpnttmp(:,3)= zlim_1;
        contr_t_top1 = [corrctrlpnttmp ones(Npts,1)]*tform';

        % second layer at z=25
        zlim_2 = zlimdefaults(2);
        corrctrlpnttmp(:,3)= zlim_2;
        contr_t_top2 = [corrctrlpnttmp ones(Npts,1)]*tform';

        % third layer at z=251-26
        zlim_3 = zlimdefaults(3);
        corrctrlpnttmp(:,3)= zlim_3;
        contr_t_bot1 = [corrctrlpnttmp ones(Npts,1)]*tform';

        % forth layer at z=251-16
        zlim_4 = zlimdefaults(4);
        corrctrlpnttmp(:,3)= zlim_4;
        contr_t_bot2 = [corrctrlpnttmp ones(Npts,1)]*tform';

        topcntrl1(:,:,tile_index) = contr_t_top1;
        topcntrl2(:,:,tile_index) = contr_t_top2;
        botcntrl1(:,:,tile_index) = contr_t_bot1;
        botcntrl2(:,:,tile_index) = contr_t_bot2;

        control(1:2*Npts,:,tile_index) = [contr_t_top1;contr_t_top2];
        control(2*Npts+1:end,:,tile_index) = [contr_t_bot1;contr_t_bot2];
        zlim_cntrl(:,tile_index) = [zlim_1 zlim_2 zlim_3 zlim_4];

        % pix stats
        layer = regpts{tile_index}.X;
        layerp1 = regpts{tile_index}.Y;
        if isempty(layer)
            continue
        end
        meddesc = round(median([layer(:,3) layerp1(:,3)],1));
        mindesc = round(min([layer(:,3) layerp1(:,3)],[],1));
        maxdesc = round(max([layer(:,3) layerp1(:,3)],[],1));
        pixstats(tile_index,:) = [tile_index  regpts{tile_index}.neigs(4) meddesc mindesc maxdesc];
        % pixstats(62,2:end) = pixstats(1162,2:end);
    end

    Rmin = min(scopeloc.loc)*1e6;
    Rmin=Rmin*.995;
    Rmax = round((max(scopeloc.loc)*1e3+scopeparams(1).imsize_um)*1e3);
    Rmax=Rmax*1.005;
    targetidx = [];
    latticeZRange = unique(scopeloc.gridix(:,3));
    if nargin<6 || isempty(theselayers) ,
        theselayers=latticeZRange(1:end-1)';
    end
    nooptim = NaN(1,tile_count);
    for t = theselayers ,
        %
        fprintf(['    Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3))) '\n']);
        ix = (scopeloc.gridix(:,3)'==t);
        idxinlayer = find(ix);
        if sum(ix)<1 ,
            fprintf(['No tiles found in layer #' num2str(t) '!!\n']);
            continue
        end
        
        % get interpolants based on paired descriptors
        [Fxt, Fyt, Fzt, Fxtp1, Fytp1, Fztp1,XYZ_tori,XYZ_tp1ori,outliers] =...
            util.getInterpolants(ix,regpts,afftile,params,curvemodel);
        
        %
        if params.debug ,
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
        
        %
        if isempty(Fxt) || size(Fxt.Points,1)<10 ,
            fprintf(['    MISSING SLICE @ Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3))) '\n']);
            continue,
        else
            fprintf(['    Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3))) ' totDesc: ' num2str(size(Fxt.Points,1)) '\n']);
        end
        
        %
        for tile_index = idxinlayer ,  % layer t
            if ~isempty(targetidx)
                if ~any(tile_index==targetidx)
                    continue
                end
            end
            %idxtm1 = tileneighbors(tile_index,6);
            idxtp1 = tileneighbors(tile_index,7);

            %
            locs = util.fcshift(curvemodel(:,:,tile_index),order,xy,dims,subcorrctrlpnttmp) ;
            corrctrlpnttmp = locs-1;

            [nooptim(tile_index),control_t_bot12,control_tp1_top12,zlim_cntrl] = optimpertile...
                (tile_index,params,tileneighbors,afftile,pixstats,zlim_cntrl,corrctrlpnttmp,Fxt,Fyt,Fzt,Fxtp1,Fytp1,Fztp1);
            if nooptim(tile_index) % nomatch with below

            elseif isempty(control_tp1_top12) % no below adjacent tile
                control(2*Npts+1:end,:,tile_index) = control_t_bot12;
            else
                control(2*Npts+1:end,:,tile_index) = control_t_bot12;
                control(1:2*Npts,:,idxtp1) = control_tp1_top12;
            end

        end
        
        %
        anchorinds = idxinlayer(~nooptim(idxinlayer));
        anchors = scopeloc.gridix(anchorinds,:);        
        queryinds = idxinlayer(nooptim(idxinlayer)>0);
        queries = scopeloc.gridix(queryinds,:);
        IDX_nn = knnsearch(anchors,queries,'K',1);
        IDX = rangesearch(anchors,queries,sqrt(2));
        
        %
        for ine = 1 : length(queryinds) ,
            if isempty(IDX{ine})
                IDX{ine} = IDX_nn(ine);
            end
        end
        
        %
        for tile_index = idxinlayer ,  % layer t
            % idxtm1 = tileneighbors(tile_index,6);
            idxtp1 = tileneighbors(tile_index,7);
            if nooptim(tile_index)
                is_query =(tile_index==queryinds) ;
                zlimdefaults = round(median(zlim_cntrl(:,anchorinds(IDX{is_query})),2))';
                
                [nooptim(tile_index),control_t_bot12,control_tp1_top12,zlim_cntrl] = ...
                    nomatchoptim(tile_index,params,tileneighbors,afftile,pixstats,zlim_cntrl,corrctrlpnttmp,...
                                 Fxt,Fyt,Fzt,Fxtp1,Fytp1,Fztp1,zlimdefaults);

                control(2*Npts+1:end,:,tile_index) = control_t_bot12;
                control(1:2*Npts,:,idxtp1) = control_tp1_top12;
            end
        end
    end
    bbox = zeros(8,3,tile_count);
    origin = zeros(tile_count,3);
    sz = zeros(tile_count,3);
    for i = 1:tile_count
        [bbox(:,:,i), origin(i,:), sz(i,:)] = util.bboxFromCorners(control(:,:,i));
    end
    
    % Update per-tile affine transforms based on control points
    fprintf('Updating per-tile affine transforms based on control points...\n')
    tform = zeros(5,5,tile_count);
    numX = size(control(:,:,1),1);
    for tile_index = 1 : tile_count ,
        X = control(:,:,tile_index)'/1000;
        x = xlim_cntrl;
        y = ylim_cntrl;
        z = zlim_cntrl(:,tile_index)';
        [xx,yy,zz] = ndgrid(x,y,z);
        r = [xx(:),yy(:),zz(:)]';
        X_aug = [X;ones(1,numX)];
        r_aug = [r;ones(1,numX)];
        %A = (X*x')/(x*x');
        A = X_aug/r_aug ;
        Aest = eye(5);
        Aest(1:3,1:3) = A(1:3,1:3)*1000;
        Aest(1:3,5) = A(1:3,4)*1000;
        Aest(5,1:3) = A(4,1:3)*1000;
        tform(:,:,tile_index) = Aest;
    end
    fprintf('Done updating per-tile affine transforms based on control points.\n')
    
    % Package everything up in a single struct for return
    vecfield = struct() ;
    vecfield.control = control;
    vecfield.xlim_cntrl = xlim_cntrl;
    vecfield.ylim_cntrl = ylim_cntrl;
    vecfield.zlim_cntrl = zlim_cntrl;
    vecfield.afftile = afftile;
    vecfield.tform = tform;
    vecfield.corrctrlpnttmp = corrctrlpnttmp;
    vecfield.ctrlpnttmp = subcorrctrlpnttmp-1;
    vecfield.bbox = bbox;
    vecfield.origin = origin;
    vecfield.sz = sz;
    vecfield.time = datestr(now);
    vecfield.theselayers = theselayers;
end
