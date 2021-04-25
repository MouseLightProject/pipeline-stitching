function [Fx, Fy, Fz, Fx_neighbor, Fy_neighbor, Fz_neighbor, layer_matches_xyz, neighbor_layer_matches_xyz, outliers] = ...
    getInterpolants(tile_index_from_tile_within_layer_index, ...
                    regpts, ...
                    affine_transform_from_tile_index, ...
                    tile_ij1s, ...
                    params, ...
                    curve_model, ...
                    do_apply_field_correction)

    % Break out parameters
    do_visualize = params.viz;
    tile_shape_ijk = params.imagesize ;
    if isfield(params, 'order') ,
        order = params.order ;
    else
        order = 1 ;
    end

    % Collect up all the matched landmarks from this layer and the z+1 layer    
    tiles_with_enough_matches_count = 0 ;
    tile_within_layer_count = length(tile_index_from_tile_within_layer_index) ;
    tile_matches_xyz_from_tile_within_layer_index = cell(1, tile_within_layer_count) ;
    neighbor_matches_xyz_from_tile_within_layer_index = cell(1, tile_within_layer_count) ;
    for tile_within_layer_index = 1 : tile_within_layer_count ,
        tile_index = tile_index_from_tile_within_layer_index(tile_within_layer_index) ;

        % Get the matches with the z+1 tile
        this_tile_regpts = regpts{tile_index} ;
        raw_tile_matches_ijk0 = this_tile_regpts.X ;
        raw_neighbor_matches_ijk0 = this_tile_regpts.Y ;
        match_count = size(raw_tile_matches_ijk0, 1) ;
        
        % If not enough matches, bail on this tile
        if match_count < 250 ,
            tile_matches_xyz_from_tile_within_layer_index{tile_within_layer_index} = zeros(0,3) ;
            neighbor_matches_xyz_from_tile_within_layer_index{tile_within_layer_index} = zeros(0,3) ;            
        else
            tiles_with_enough_matches_count = tiles_with_enough_matches_count + 1 ;

            % Get the tile index of the z+1 neighbor tile
            tile_index_from_neighbor_index = this_tile_regpts.neigs ;
            self_tile_index = tile_index_from_neighbor_index(1) ;
            if self_tile_index ~= tile_index ,
                error('Something has gone horribly wrong: the self tile index (%d) does not equal the tile_index (%d)!\n', self_tile_index, tile_index) ;
            end
            neighbor_tile_index = tile_index_from_neighbor_index(4) ;

            % Apply field correction if called for
            if do_apply_field_correction ,
                this_curve_model = curve_model(:,:,tile_index) ;
                raw_tile_matches_ijk1 = raw_tile_matches_ijk0 + 1 ;
                tile_matches_ijk1 = util.fcshift(this_curve_model, order, tile_ij1s, tile_shape_ijk, raw_tile_matches_ijk1) ;
                tile_matches_ijk0 = tile_matches_ijk1 - 1 ;
                neighbor_curve_model = curve_model(:,:,neighbor_tile_index) ;
                raw_neighbor_matches_ijk1 = raw_neighbor_matches_ijk0 + 1 ;
                neighbor_matches_ijk1 = util.fcshift(neighbor_curve_model, order, tile_ij1s, tile_shape_ijk, raw_neighbor_matches_ijk1) ;
                neighbor_matches_ijk0 = neighbor_matches_ijk1 - 1 ;
            else
                tile_matches_ijk0 = raw_tile_matches_ijk0 ;
                neighbor_matches_ijk0 = raw_neightbor_matches_ijk0 ;            
            end

            % Apply per-tile affine transforms
            tile_affine_transform = affine_transform_from_tile_index(:,:,tile_index) ;
            tile_matches_xyz = add_ones_column(tile_matches_ijk0) * tile_affine_transform' ;
            neighbor_affine_transform = affine_transform_from_tile_index(:,:,neighbor_tile_index) ;
            neighbor_matches_xyz = add_ones_column(neighbor_matches_ijk0) * neighbor_affine_transform' ;

            tile_matches_xyz_from_tile_within_layer_index{tile_within_layer_index} = tile_matches_xyz ;
            neighbor_matches_xyz_from_tile_within_layer_index{tile_within_layer_index} = neighbor_matches_xyz ;
        end
    end

    % Collect all the matches in the layer, and the z+1 layer
    layer_matches_xyz = cat(1,tile_matches_xyz_from_tile_within_layer_index{:}) ;
    neighbor_layer_matches_xyz = cat(1,neighbor_matches_xyz_from_tile_within_layer_index{:}) ;
    
    % If there are no tiles with enough matches in the whole layer, just exit now, returning everything
    % empty
    layer_match_count = size(layer_matches_xyz, 1) ;
    if layer_match_count == 0 , 
        Fx = [] ;
        Fy = [] ;
        Fz = [] ;
        Fx_neighbor = [] ;
        Fy_neighbor = [] ;
        Fz_neighbor = [] ;
        layer_matches_xyz = [] ;
        neighbor_layer_matches_xyz = [] ;
        outliers = [] ;            
        return
    end
    
    %
    % eliminate outliers
    %
    
    % build kd-tree based on xy location
    K = min(20, round(sqrt(layer_match_count))) ;
    layer_matches_xy = layer_matches_xyz(:,1:2) ;
    IDX = knnsearch(layer_matches_xy, layer_matches_xy, 'K', K) ;
%     if 1
        diffXY = layer_matches_xyz(:,1:2)-neighbor_layer_matches_xyz(:,1:2);
        normdiffXY = normr(diffXY);
        % interpolate vector from nearest K samples
        bins = [linspace(-1,1,21)-.05 1.05];
        bins2 = [linspace(0.1,1,10)-.05 1.05];
        st = zeros(size(diffXY,1),4);
        for idx = 1:size(diffXY,1)
            dists = [0;sqrt(sum((ones(size(IDX,2)-1,1)*layer_matches_xyz(idx,1:2)-layer_matches_xyz(IDX(idx,2:end),1:2)).^2,2))]; % can be used as weighting
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
        outliers1 = st(:,2)<.8;
        outliers2 = st(:,3)<.5 & st(:,4)<.5;
        outliers = outliers1|outliers2;
        inliers = ~outliers;
%     else
%         diffZ = layer_matches_xyz(:,3)-neighbor_layer_matches_xyz(:,3);diffZ=diffZ(IDX);
%         % compare first to rest
%         inliers = abs(diffZ(:,1))<=2e3 | abs(diffZ(:,1)) <= 3*abs(median(diffZ(:,2:end),2));
%     end
    
    if do_visualize
        figure(35), cla
        hold on
        plot3(layer_matches_xyz(:,1),layer_matches_xyz(:,2),layer_matches_xyz(:,3),'r.') % layer t
        plot3(neighbor_layer_matches_xyz(:,1),neighbor_layer_matches_xyz(:,2),neighbor_layer_matches_xyz(:,3),'k.') % layer tp1
        %text(XYZ_t(outliers,1),XYZ_t(outliers,2),num2str(st(outliers,2:4)))
        myplot3(layer_matches_xyz(outliers,:),'md') % layer t
        myplot3(neighbor_layer_matches_xyz(outliers,:),'gd') % layer t
        %plot3(XYZ_t(IDX(idx,1),1),XYZ_t(IDX(idx,1),2),XYZ_t(IDX(idx,1),3),'r*') % layer t
        %plot3(XYZ_t(IDX(idx,2:end),1),XYZ_t(IDX(idx,2:end),2),XYZ_t(IDX(idx,2:end),3),'bo') % layer t
        set(gca,'Ydir','reverse')
        drawnow
    end
    
    inlier_layer_matches_xyz = layer_matches_xyz(inliers,:);
    inlier_neighbor_layer_matches_xyz = neighbor_layer_matches_xyz(inliers,:);
    
    if do_visualize
        figure(33), cla
        hold on
        plot3(inlier_layer_matches_xyz(:,1),inlier_layer_matches_xyz(:,2),inlier_layer_matches_xyz(:,3),'k.') % layer t
        plot3(inlier_neighbor_layer_matches_xyz(:,1),inlier_neighbor_layer_matches_xyz(:,2),inlier_neighbor_layer_matches_xyz(:,3),'r.') % layer tp1
        myplot3([tile_matches_ijk0 ones(size(tile_matches_ijk0,1),1)]*affine_transform_from_tile_index(:,:,tile_index)','g.') % layer t
        myplot3([neighbor_matches_ijk0 ones(size(tile_matches_ijk0,1),1)]*affine_transform_from_tile_index(:,:,neighbor_tile_index)','m.') % layer tp1
        set(gca,'Zdir','reverse');
        legend('layer t','layer t+1')
        set(gca,'Ydir','reverse')
        %         plot3(XYZ_t(inliers,1),XYZ_t(inliers,2),XYZ_t(inliers,3),'go') % layer t
        %         plot3(XYZ_tp1(inliers,1),XYZ_tp1(inliers,2),XYZ_tp1(inliers,3),'mo') % layer tp1
        %         plot3(XYZ_t(16248,1),XYZ_t(16248,2),XYZ_t(16248,3),'g*') % layer t
        %         plot3(XYZ_tp1(16248,1),XYZ_tp1(16248,2),XYZ_tp1(16248,3),'m*') % layer tp1
    end
    
    vecdif = inlier_layer_matches_xyz - inlier_neighbor_layer_matches_xyz ;
    rt = params.expensionratio/(1+params.expensionratio) ;
    
    % Get an interpolation function for x, y, and z for this layer
    Fx = scatteredInterpolant(inlier_layer_matches_xyz, -vecdif(:,1)*(1-rt), 'linear', 'nearest') ;
    Fy = scatteredInterpolant(inlier_layer_matches_xyz, -vecdif(:,2)*(1-rt), 'linear', 'nearest') ;
    Fz = scatteredInterpolant(inlier_layer_matches_xyz, -vecdif(:,3)*(1-rt), 'linear', 'nearest') ;
    
    % Get an interpolation function for x, y, and z for the z+1 layer
    Fx_neighbor = scatteredInterpolant(inlier_neighbor_layer_matches_xyz, vecdif(:,1)*rt, 'linear', 'nearest') ;
    Fy_neighbor = scatteredInterpolant(inlier_neighbor_layer_matches_xyz, vecdif(:,2)*rt, 'linear', 'nearest') ;
    Fz_neighbor = scatteredInterpolant(inlier_neighbor_layer_matches_xyz, vecdif(:,3)*rt, 'linear', 'nearest') ;
end
