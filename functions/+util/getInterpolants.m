function [Fx_layer, Fy_layer, Fz_layer, Fx_next_layer, Fy_next_layer, Fz_next_layer, layer_xyz_from_match_index, next_layer_xyz_from_match_index, is_outlier_from_match_index] = ...
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
    expansion_ratio = params.expensionratio ;
    if isfield(params, 'order') ,
        order = params.order ;
    else
        order = 1 ;
    end

    % Collect up all the matched landmarks from this layer and the z+1 layer
    [layer_xyz_from_match_index, next_layer_xyz_from_match_index] = ...
        collect_layer_matches(tile_index_from_tile_within_layer_index, ...
                              regpts, ...
                              curve_model, ...
                              order, ...
                              affine_transform_from_tile_index, ...
                              tile_shape_ijk, ...
                              tile_ij1s, ...
                              do_apply_field_correction) ;
    
    % If there are no tiles with enough matches in the whole layer, just exit now, returning everything
    % empty
    match_count = size(layer_xyz_from_match_index, 1) ;
    if match_count == 0 , 
        Fx_layer = [] ;
        Fy_layer = [] ;
        Fz_layer = [] ;
        Fx_next_layer = [] ;
        Fy_next_layer = [] ;
        Fz_next_layer = [] ;
        layer_xyz_from_match_index = [] ;
        next_layer_xyz_from_match_index = [] ;
        is_outlier_from_match_index = [] ;            
        return
    end
    
    % Identify outliers
    % Do this by computing a set of features for each match, all of which we'd like
    % to be high (close to one).  If enough of them are much less than one, we
    % consider the match to be an outlier.
    K = min(20, round(sqrt(match_count))) ;
    layer_xy_from_match_index = layer_xyz_from_match_index(:,1:2) ;
    next_layer_xy_from_match_index = next_layer_xyz_from_match_index(:,1:2) ;
    IDX = knnsearch(layer_xy_from_match_index, layer_xy_from_match_index, 'K', K) ;
    
    % Do a sanity check, b/c the code used to use IDX(match_index,1) in several
    % places where it seems match_index would do...
    if ~isequal(IDX(:,1), (1:match_count)') ,
        error('Something has gone horribly wrong: IDX(:,1) ~= (1:match_count)''') ;
    end            
    
    dxy_from_match_index = layer_xy_from_match_index - next_layer_xy_from_match_index ;  
      % All the "dxy" variables are differences in the matched points between the two
      % layers, relative to the next layer
    dxy_hat_from_match_index = normr(dxy_from_match_index) ;
    % interpolate vector from nearest K samples
    similarity_bin_edges = -1.05 : 0.1 : 1.05 ;
    similarity_bin_centers = (similarity_bin_edges(1:end-1) + similarity_bin_edges(2:end))/2 ;
    positive_bin_edges = -0.05 : 0.1 : 1.05 ;
    positive_bin_centers = (positive_bin_edges(1:end-1) + positive_bin_edges(2:end))/2 ;
    modal_similarity_from_match_index = zeros(match_count,1) ;
    modal_boltzman_similarity_from_match_index = zeros(match_count,1) ;
    approximate_mean_boltzman_similarity_from_match_index = zeros(match_count,1) ;    
    for match_index = 1 : match_count ,
        % Get a bunch of info about this match, and nearby matches
        layer_xy = layer_xy_from_match_index(match_index,:) ;
        match_index_from_nearby_match_index = IDX(match_index,:) ;  % includes this_match_xy in row 1
        layer_xy_from_nearby_match_index = layer_xy_from_match_index(match_index_from_nearby_match_index, :) ;
        dxy = dxy_from_match_index(match_index,:) ;
        dxy_from_nearby_match_index = dxy_from_match_index(match_index_from_nearby_match_index,:) ;
        dxy_hat = dxy_hat_from_match_index(match_index,:) ;
        dxy_hat_from_nearby_match_index = dxy_hat_from_match_index(match_index_from_nearby_match_index,:) ;

        % Make a filter for which ones are at an acceptable distance
        xy_distance_from_nearby_match_index = vecnorm(layer_xy_from_nearby_match_index - layer_xy, 2, 2) ;  % 1 x nearby_match_count
            % distance in xy from each nearby match to the current match, all within the
            % current layer
        is_at_good_distance = ( 0<xy_distance_from_nearby_match_index & xy_distance_from_nearby_match_index<1e6 ) ;  % 1e6 is in nm
        is_at_good_distance(1) = false ;  % make double-sure match_index is not included
                
        % Find the approximate mode in the distribution of
        % cosine similarity between this match and nearby matches.
        similarity_from_nearby_match_index = dxy_hat * dxy_hat_from_nearby_match_index' ;  
            % cosine similarity in x-y plane from each nearby match to the current match,
            % all within the current layer
        modal_similarity = ...
            mode_after_binning(similarity_from_nearby_match_index(is_at_good_distance), ...
                               similarity_bin_edges, ...
                               similarity_bin_centers) ;  % "mode similarity", akin to "mean similarity" or "median similarity"             
        modal_similarity_from_match_index(match_index) = modal_similarity ;

        % Compute how well this shift matches nearby shifts, in a different way
        relative_dxy_from_nearby_match_index = dxy_from_nearby_match_index - dxy ;
        relative_dxy_length_from_nearby_match_index = vecnorm(relative_dxy_from_nearby_match_index, 2, 2) ; 
             % 2nd arg specifies euclidian norm, 1 x nearby_match_count, 1st element always 0
        boltzman_similarity_from_nearby_match_index = exp(-relative_dxy_length_from_nearby_match_index/norm(dxy)) ;  % 1 x nearby_match_count, 1st element always 1
            % exp(-x) converts a distance on [0,inf) to a similarity measure on (0,1]
        modal_boltzman_similarity = ...
            mode_after_binning_favoring_high(boltzman_similarity_from_nearby_match_index(is_at_good_distance), ...
                                             positive_bin_edges, ...
                                             positive_bin_centers) ;  % "mode similarity", akin to "mean similarity" or "median similarity"             
        modal_boltzman_similarity_from_match_index(match_index) = modal_boltzman_similarity ;

        boltzman_similarity_count_from_bin_index = ...
            histc(boltzman_similarity_from_nearby_match_index(is_at_good_distance), positive_bin_edges) ;  %#ok<HISTC>  % Max-likely
        approximate_mean_boltzman_similarity = (positive_bin_edges-.05)*boltzman_similarity_count_from_bin_index/sum(boltzman_similarity_count_from_bin_index) ;  
          % seems like an approximation of the mean of boltzman_similarity_from_nearby_match_index(is_at_good_distance)
          % but (positive_bin_edges-.05) should be (positive_bin_edges+.05).
          % Or just do
          % mean(boltzman_similarity_from_nearby_match_index(is_at_good_distance))...
        approximate_mean_boltzman_similarity_from_match_index(match_index) = approximate_mean_boltzman_similarity ;
    end
    is_type_1_outlier_from_match_index = modal_similarity_from_match_index < .8 ;
    is_type_2_outlier_from_match_index = modal_boltzman_similarity_from_match_index<.5 & approximate_mean_boltzman_similarity_from_match_index<.5 ;
        % Not like type I and type II errors, just different reasons to exclude
    is_outlier_from_match_index = is_type_1_outlier_from_match_index | is_type_2_outlier_from_match_index ;
    
    % Visualize the matches (inliers and outliers)
    if do_visualize
        figure(35), cla
        hold on
        plot3(layer_xyz_from_match_index(:,1),layer_xyz_from_match_index(:,2),layer_xyz_from_match_index(:,3),'r.') % layer t
        plot3(next_layer_xyz_from_match_index(:,1),next_layer_xyz_from_match_index(:,2),next_layer_xyz_from_match_index(:,3),'k.') % layer tp1
        %text(XYZ_t(outliers,1),XYZ_t(outliers,2),num2str(st(outliers,2:4)))
        myplot3(layer_xyz_from_match_index(is_outlier_from_match_index,:),'md') % layer t
        myplot3(next_layer_xyz_from_match_index(is_outlier_from_match_index,:),'gd') % layer t
        %plot3(XYZ_t(IDX(idx,1),1),XYZ_t(IDX(idx,1),2),XYZ_t(IDX(idx,1),3),'r*') % layer t
        %plot3(XYZ_t(IDX(idx,2:end),1),XYZ_t(IDX(idx,2:end),2),XYZ_t(IDX(idx,2:end),3),'bo') % layer t
        set(gca,'Ydir','reverse')
        drawnow
    end
    
    % Discard outliers
    is_inlier_from_match_index = ~is_outlier_from_match_index ;
    layer_xyz_from_inlier_match_index = layer_xyz_from_match_index(is_inlier_from_match_index,:) ;
    next_layer_xyz_from_inlier_match_index = next_layer_xyz_from_match_index(is_inlier_from_match_index,:) ;
    
    % Visualize the inliers
    if do_visualize
        figure(33), cla
        hold on
        plot3(layer_xyz_from_inlier_match_index(:,1),layer_xyz_from_inlier_match_index(:,2),layer_xyz_from_inlier_match_index(:,3),'k.') % layer t
        plot3(next_layer_xyz_from_inlier_match_index(:,1),next_layer_xyz_from_inlier_match_index(:,2),next_layer_xyz_from_inlier_match_index(:,3),'r.') % layer tp1
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
    
    % Compute the shifts implied by the matches, and the how much each tile will
    % account for the shift
    dxyz = layer_xyz_from_inlier_match_index - next_layer_xyz_from_inlier_match_index ;
    neighbor_weight = expansion_ratio/(1+expansion_ratio) ;  % expansion_ratio is usually one, so this is usually 1/2
    tile_weight = 1 - neighbor_weight ;
    tile_dxyz = -tile_weight * dxyz ;
    neighbor_dxyz = neighbor_weight * dxyz ;
    
    % Get an interpolation function for x, y, and z for this layer
    Fx_layer = scatteredInterpolant(layer_xyz_from_inlier_match_index, tile_dxyz(:,1), 'linear', 'nearest') ;
    Fy_layer = scatteredInterpolant(layer_xyz_from_inlier_match_index, tile_dxyz(:,2), 'linear', 'nearest') ;
    Fz_layer = scatteredInterpolant(layer_xyz_from_inlier_match_index, tile_dxyz(:,3), 'linear', 'nearest') ;
    
    % Get an interpolation function for x, y, and z for the z+1 layer
    Fx_next_layer = scatteredInterpolant(next_layer_xyz_from_inlier_match_index, neighbor_dxyz(:,1), 'linear', 'nearest') ;
    Fy_next_layer = scatteredInterpolant(next_layer_xyz_from_inlier_match_index, neighbor_dxyz(:,2), 'linear', 'nearest') ;
    Fz_next_layer = scatteredInterpolant(next_layer_xyz_from_inlier_match_index, neighbor_dxyz(:,3), 'linear', 'nearest') ;
end



function [layer_matches_xyz, neighbor_layer_matches_xyz] = ...
        collect_layer_matches(tile_index_from_tile_within_layer_index, ...
                              regpts, ...
                              curve_model, ...
                              order, ...
                              affine_transform_from_tile_index, ...
                              tile_shape_ijk, ...
                              tile_ij1s, ...
                              do_apply_field_correction)
    tile_within_layer_count = length(tile_index_from_tile_within_layer_index) ;
    tile_matches_xyz_from_tile_within_layer_index = cell(1, tile_within_layer_count) ;
    neighbor_matches_xyz_from_tile_within_layer_index = cell(1, tile_within_layer_count) ;
    for tile_within_layer_index = 1 : tile_within_layer_count ,
        tile_index = tile_index_from_tile_within_layer_index(tile_within_layer_index) ;

        % Get the matches with the z+1 tile
        this_tile_regpts = regpts{tile_index} ;
        raw_tile_matches_ijk0 = this_tile_regpts.X ;
        raw_neighbor_matches_ijk0 = this_tile_regpts.Y ;
        tile_match_count = size(raw_tile_matches_ijk0, 1) ;
        
        % If not enough matches, bail on this tile
        if tile_match_count < 250 ,
            tile_matches_xyz_from_tile_within_layer_index{tile_within_layer_index} = zeros(0,3) ;
            neighbor_matches_xyz_from_tile_within_layer_index{tile_within_layer_index} = zeros(0,3) ;            
        else
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
end



function result = mode_after_binning(value_from_index, bin_edges, bin_centers)
    % Find the approximate mode of a set of data by making a histogram and returning
    % the bin center of the bin with the highest count.
    % Typically, bin_centers will be pre-computed according to
    %   bin_centers = ( bin_edges(1:end-1) + bin_edges(2:end) ) / 2
    % In any case, bin_edges should be one element longer than bin_centers.
    counts_from_bin_index = histc(value_from_index, bin_edges) ;  %#ok<HISTC>
    counts_from_bin_index = counts_from_bin_index(1:end-1) ;  % drop the last "bin", since it's the zero-width one
    [~,modal_bin_index] = max(counts_from_bin_index) ;  % defaults to the lowest argmax if there are ties
    result = bin_centers(modal_bin_index) ;
end



function result = mode_after_binning_favoring_high(value_from_index, bin_edges, bin_centers)
    % Find the approximate mode of a set of data by making a histogram and returning
    % the bin center of the bin with the highest count.
    % Typically, bin_centers will be pre-computed according to
    %   bin_centers = ( bin_edges(1:end-1) + bin_edges(2:end) ) / 2
    % In any case, bin_edges should be one element longer than bin_centers.
    % If bins are tied for the mode, this version returns the center of the
    % highest-value bin.
    counts_from_bin_index = histc(value_from_index, bin_edges) ;  %#ok<HISTC>
    counts_from_bin_index = counts_from_bin_index(1:end-1) ;  % drop the last "bin", since it's the zero-width one
    max_counts = max(counts_from_bin_index) ;
    is_mode_from_bin_index = (counts_from_bin_index==max_counts) ;
    modal_bin_index = find(is_mode_from_bin_index, 1, 'last') ;  % pick the last one if there are multiples
    result = bin_centers(modal_bin_index) ;
end
