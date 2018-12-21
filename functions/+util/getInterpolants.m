function [Fxt, Fyt, Fzt, Fxtp1, Fytp1, Fztp1,XYZ_tori,XYZ_tp1ori,outliers] = ...
    getInterpolants(ix,regpts,afftile,params,curvemodel)

viz = params.viz;
dims = params.imagesize;
cnt = 0;
Npts = 0;
if isfield(params,'order')
    order = params.order;
else
    order = 1;
end
if isfield(params,'applyFC')
    applyFC = params.applyFC;
else
    applyFC = 0;
end

xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];

[poslayer_t,poslayer_tp1,Fxt, Fyt, Fzt, Fxtp1, Fytp1, Fztp1,XYZ_tori,XYZ_tp1ori,outliers] = deal([]);
for id_ix = find(ix)
    %%
    if isempty(regpts{id_ix}.X) | size(regpts{id_ix}.X,1)<250
        continue
    end
    %%
    neigs = regpts{id_ix}.neigs;
    layer = regpts{id_ix}.X;
    layerp1 = regpts{id_ix}.Y;

    % apply FC
    if applyFC
        [layer] = util.fcshift(curvemodel(:,:,neigs(1)),order,xy,dims,layer+1);
        layer = layer-1;
        [layerp1] = util.fcshift(curvemodel(:,:,neigs(4)),order,xy,dims,layerp1+1);
        layerp1 = layerp1-1;
    end
    %%
    % apply tforms
    layer_t = [layer ones(size(layer,1),1)]*afftile(:,:,neigs(1))';
    layer_tp1 = [layerp1 ones(size(layerp1,1),1)]*afftile(:,:,neigs(4))';
    
    cnt = cnt+1;
    Npts = Npts+size(layer_t,1);
    poslayer_t{cnt} = layer_t;
    poslayer_tp1{cnt} = layer_tp1;
end
if isempty(poslayer_t),return,end
%%
XYZ_t = cat(1,poslayer_t{:});
XYZ_tp1 = cat(1,poslayer_tp1{:});
% eliminate outliers
% build kd-tree based on xy locatio
K = min(20,round(sqrt(size(XYZ_t,1))));
IDX = knnsearch(XYZ_t(:,1:2),XYZ_t(:,1:2),'K',K);
%
if 1
    diffXY = XYZ_t(:,1:2)-XYZ_tp1(:,1:2);
    normdiffXY = normr(diffXY);
    % interpolate vector from nearest K samples
    bins = [linspace(-1,1,21)-.05 1.05];
    bins2 = [linspace(0.1,1,10)-.05 1.05];
    st = zeros(size(diffXY,1),4);
    for idx = 1:size(diffXY,1)
        %%
        
        dists = [0;sqrt(sum((ones(size(IDX,2)-1,1)*XYZ_t(idx,1:2)-XYZ_t(IDX(idx,2:end),1:2)).^2,2))]; % can be used as weighting
        %         [bandwidth,density,xmesh,cdf]=kde(dists);
        
        if 1
            innprod = [1 normdiffXY(idx,:)*normdiffXY(IDX(idx,2:end),:)'];
        else
            weights = exp(-dists/std(dists(dists>0 & dists<1e6)));weights = weights/max(weights);
            innprod = [normdiffXY(idx,:)*normdiffXY(IDX(idx,:),:)'*diag(weights)];
        end
        % majority binning
        % theta
        [~,idxmaxtheta] = max(histc(innprod(dists>0 & dists<1e6),bins));
        st(idx,2) = bins(idxmaxtheta)+.05;
        
        dV = ones(size(IDX(idx,:),2),1)*diffXY(IDX(idx,1),:)-diffXY(IDX(idx,1:end),:);
        mags = sqrt(sum(dV.^2,2));
        mags = exp(-mags/norm(diffXY(IDX(idx,1),:)));
        mags = mags/max(mags);
        xx=flipud(histc(mags(dists>0 & dists<1e6),bins2));
        [tr,idxmaxdist] = max(xx); % Max-likely
        idxmaxdist = length(xx)+1-idxmaxdist;
        [aa,bb] = histc(mags(dists>0 & dists<1e6),bins2); % Max-likely
        
        st(idx,3) = bins2(idxmaxdist)+.05;
        st(idx,4) = (bins2-.05)*aa(:)/sum(aa);
    end
    %%
    outliers1 = st(:,2)<.8;
    outliers2 = st(:,3)<.5 & st(:,4)<.5;
    outliers = outliers1|outliers2;
%     outliers = outliers1;
    inliers = ~outliers;
else
    %%
    diffZ = XYZ_t(:,3)-XYZ_tp1(:,3);diffZ=diffZ(IDX);
    % compare first to rest
    inliers = abs(diffZ(:,1))<=2e3 | abs(diffZ(:,1)) <= 3*abs(median(diffZ(:,2:end),2));
end
    XYZ_tori = XYZ_t;
    XYZ_tp1ori = XYZ_tp1;
if viz
    %%
    figure(35), cla
    hold on
    plot3(XYZ_t(:,1),XYZ_t(:,2),XYZ_t(:,3),'r.') % layer t
    plot3(XYZ_tp1(:,1),XYZ_tp1(:,2),XYZ_tp1(:,3),'k.') % layer tp1
    %text(XYZ_t(outliers,1),XYZ_t(outliers,2),num2str(st(outliers,2:4)))
    myplot3(XYZ_t(outliers,:),'md') % layer t
    myplot3(XYZ_tp1(outliers,:),'gd') % layer t
    %plot3(XYZ_t(IDX(idx,1),1),XYZ_t(IDX(idx,1),2),XYZ_t(IDX(idx,1),3),'r*') % layer t
    %plot3(XYZ_t(IDX(idx,2:end),1),XYZ_t(IDX(idx,2:end),2),XYZ_t(IDX(idx,2:end),3),'bo') % layer t
    set(gca,'Ydir','reverse')
    drawnow
    
end
%%
XYZ_t = XYZ_t(inliers,:);
XYZ_tp1 = XYZ_tp1(inliers,:);
%%
if viz
    %%
    figure(33), cla
    hold on
    plot3(XYZ_t(:,1),XYZ_t(:,2),XYZ_t(:,3),'k.') % layer t
    plot3(XYZ_tp1(:,1),XYZ_tp1(:,2),XYZ_tp1(:,3),'r.') % layer tp1
    myplot3([layer ones(size(layer,1),1)]*afftile(:,:,neigs(1))','g.') % layer t
    myplot3([layerp1 ones(size(layer,1),1)]*afftile(:,:,neigs(4))','m.') % layer tp1
    set(gca,'Zdir','reverse');
    legend('layer t','layer t+1')
    set(gca,'Ydir','reverse')
    %         plot3(XYZ_t(inliers,1),XYZ_t(inliers,2),XYZ_t(inliers,3),'go') % layer t
    %         plot3(XYZ_tp1(inliers,1),XYZ_tp1(inliers,2),XYZ_tp1(inliers,3),'mo') % layer tp1
    %         plot3(XYZ_t(16248,1),XYZ_t(16248,2),XYZ_t(16248,3),'g*') % layer t
    %         plot3(XYZ_tp1(16248,1),XYZ_tp1(16248,2),XYZ_tp1(16248,3),'m*') % layer tp1
end
%%
vecdif = XYZ_t - XYZ_tp1;
rt = params.expensionratio/(1+params.expensionratio);
% layer t
Fxt = scatteredInterpolant(XYZ_t,-vecdif(:,1)*(1-rt),'linear','nearest');
Fyt = scatteredInterpolant(XYZ_t,-vecdif(:,2)*(1-rt),'linear','nearest');
Fzt = scatteredInterpolant(XYZ_t,-vecdif(:,3)*(1-rt),'linear','nearest');
% layer t+1
Fxtp1 = scatteredInterpolant(XYZ_tp1,vecdif(:,1)*rt,'linear','nearest');
Fytp1 = scatteredInterpolant(XYZ_tp1,vecdif(:,2)*rt,'linear','nearest');
Fztp1 = scatteredInterpolant(XYZ_tp1,vecdif(:,3)*rt,'linear','nearest');