function [paireddescriptor, curvemodel, unreliable, neigbors_used] = ...
        curvatureOutlierElimination(paireddescriptor, curvemodel, scopeloc, params, model)
    % compares curvature estimation wrto median model to detect outliers
    dims = params.imagesize;

    Ntile = size(curvemodel,3);
    reliable = squeeze(all(curvemodel,2) & all(isfinite(curvemodel),2))';

    % median models
    medmodels = zeros(2,3);
    for idim = 1:2
        tmp = squeeze(curvemodel(idim,:,:));
        medmodels(idim,:) = getmod(tmp(:,all(tmp))');
    end
    % compare model under extremes. 

    % baselines
    xrange = 1:dims(1);
    yrange = 1:dims(2);
    xbase = feval(model,medmodels(1,:).*[1 1 0],yrange); % 0 shift, just the curve
    ybase = feval(model,medmodels(2,:).*[1 1 0],xrange); % 0 shift

    outlier = zeros(Ntile,2); % based on initial curve fit using 75% inlier criteria
    for itile = 1:Ntile
        xest = feval(model,curvemodel(1,:,itile).*[1 1 0],yrange);
        yest = feval(model,curvemodel(2,:,itile).*[1 1 0],xrange);
        outlier(itile,1) = mean(abs(xbase-xest));
        outlier(itile,2) = mean(abs(ybase-yest));
    end
    unreliable = ~reliable(:,1:2) | outlier > 1 | isnan(outlier);


%     %
%     % outlier rejection per tile
%     %
%     inliers = find(~any(unreliable,2));
%     % for every tiles estimate an affine
%     anchors = scopeloc.gridix(inliers,1:3);
%     queries = scopeloc.gridix(:,1:3);
%     IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]
% 
%     % fill missing
%     for ineig = 1:Ntile
%         ianch = inliers(IDX(ineig));
%         paireddescriptor{ineig}.onx = paireddescriptor{ianch}.onx;
%         paireddescriptor{ineig}.ony = paireddescriptor{ianch}.ony;
%         paireddescriptor{ineig}.count = [size(paireddescriptor{ineig}.onx.X,1) size(paireddescriptor{ineig}.ony.X,1)];
%         curvemodel(:,:,ineig) = curvemodel(:,:,ianch);
%     end


    %
    % outlier rejection per direction
    %
    % util.debug.vizCurveStats() ;

    neigbors_used = zeros(Ntile,2) ;
    
    % iterate on x direction
    queries = scopeloc.gridix(:,1:3);
    inliers = find(~unreliable(:,1));
    if isempty(inliers) ,
        error('During field correction, unable to calculate a trustworthy field correction between *any* tile and it''s x+1 tile.  Giving up.') ;
    end
    anchors = scopeloc.gridix(inliers,1:3);
    IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]
    neigbors_used(:,1) = inliers(IDX);

    % fill missing on x direction
    for itile = 1:Ntile
        ianch = inliers(IDX(itile));
        paireddescriptor{itile}.onx = paireddescriptor{ianch}.onx;
        curvemodel(1,:,itile) = curvemodel(1,:,ianch);
    end

    % iterate on y direction
    inliers = find(~unreliable(:,2));
    if isempty(inliers) ,
        error('During field correction, unable to calculate a trustworthy field correction between *any* tile and it''s y+1 tile.  Giving up.') ;
    end
    anchors = scopeloc.gridix(inliers,1:3);
    IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]
    neigbors_used(:,2) = inliers(IDX);

    % fill missing on y direction
    for itile = 1:Ntile
        ianch = inliers(IDX(itile));
        paireddescriptor{itile}.ony = paireddescriptor{ianch}.ony;
        curvemodel(2,:,itile) = curvemodel(2,:,ianch);
    end

    for itile = 1:Ntile
        paireddescriptor{itile}.count = [size(paireddescriptor{itile}.onx.X,1) size(paireddescriptor{itile}.ony.X,1)];
    end
end
