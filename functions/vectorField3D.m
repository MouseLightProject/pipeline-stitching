function vecfield = vectorField3D(params, scopeloc, regpts, scopeparams, curvemodel, tile_k_from_run_layer_index)
    fprintf('Calculating vector fields...\n') ;
    tile_count = length(scopeloc.filepath) ;
    cpg_k_count = params.Nlayer ;  % Number of k/z values in the per-tile control point grid (traditionally 4)
    cpg_i_count = params.Ndivs+1 ;  % Number of i/x values in the per-tile control point grid (traditionally 5)
    cpg_j_count = params.Ndivs+1 ;  % Number of j/y values in the per-tile control point grid (traditionally 5)    
    cpg_ij_count = cpg_i_count * cpg_j_count ;  % Number of points in each k/z layer in the per-tile control point grid (traditionally 25)
    tile_shape_ijk = params.imagesize ;  % tile shape, in xyz order (traditionally [1024 1536 251])
    order = params.order ;  % order of the field curvature model, I think
    default_cpg_k0_values = [2 20 tile_shape_ijk(3)-20-1 tile_shape_ijk(3)-2-1] ;  
      % The z values in the control point gird (traditionally [2 20 230 248])
      % If these are symmetric, it implies that this is using zero-based indexing.
    %params.zlimdefaults = zlimdefaults;
    do_apply_field_correction = params.applyFC ;  % whether or not to apply the field correction, I think

    % Build the neighbor index
    tile_ijk_from_tile_index = scopeloc.gridix(:,1:3) ;
    tileneighbors = buildNeighbor(tile_ijk_from_tile_index) ;  % tile_count x 7, each row in [self -x -y +x +y -z +z] format
    has_k_plus_1_tile_from_tile_index = ~isnan(tileneighbors(:,7)) ;
    
    % Determine the x-y grid used for the control points of the barymetric transform
    % in the per-tile space.
    % The x-y grid a crop of the full tile x-y range, to minimize overlap between
    % tiles in the rendered space.
    [st, ed] = util.getcontolpixlocations(scopeloc, params, scopeparams) ;  % st 3x1, the first x/y/z lim in each dimension; ed 3x1, the last x/y/z lim in each dimension
    [cpg_ij0s, cpg_i0_values, cpg_j0_values] = ...
        util.getCtrlPts(tile_shape_ijk(1), tile_shape_ijk(2), params, st, ed) ;  % "cpg" is the "control point grid"
      % cpg_i0_values is traditionally 1 x 5, gives the i/x values used in the control
      % point grid.  Uses zero-based indexing
      % cpg_j0_values is traditionally 1 x 5, gives the j/y values used in the control
      % point grid.  Uses zero-based indexing
      % cpg_ij0s is traditionally 25 x 2, gives all the combinations of cpg_i0_values
      % and cpg_j0_values.  I.e. it's a raster-scan of the xy grid in the tile
    cpg_ij1s = cpg_ij0s + 1 ;  % traditionally 25 x 2, the ij/xy coords of the control point grid in each tile, in each z plane
    
    % Get the linear transform used for each tile
    % Need to flip the x/y dimensions in the per-tile linear transform to match
    % imaging
    raw_linear_transform_from_tile_index = reshape([scopeparams(:).affineglFC],3,3,[]);
    mflip = -eye(3) ;
    linear_transform_from_tile_index = zeros(size(raw_linear_transform_from_tile_index)) ;
    for ii = 1 : size(raw_linear_transform_from_tile_index,3) ,
        linear_transform_from_tile_index(:,:,ii) = raw_linear_transform_from_tile_index(:,:,ii) * mflip ;
    end

    % For efficiency reasons, precompute the grid of i,j values in each plane of the
    % tile (for all pixels, not on the control point grid)
    xlocs = 1:tile_shape_ijk(1) ;
    ylocs = 1:tile_shape_ijk(2) ;
    [xy2, xy1] = ndgrid(ylocs(:), xlocs(:)) ;
    tile_ij1s = [xy1(:), xy2(:)] ;  % tile_shape_ijk(1)*tile_shape_ijk(2) x 2, the ij1/xy coords of each pixel in a z plane of the tile stack

    % Extract the root of the raw tiles folder 
    filep = strsplit(scopeloc.filepath{1},'/') ;
    vecfield_root = fullfile('/',filep{1:end-4}) ;  % e.g. '/groups/mousebrainmicro/mousebrainmicro/data/2020-09-15/Tiling'

    % Get the relative path to each tile from scopeloc
    vecfield_path = cell(1, tile_count) ;
    for tile_index = 1 : tile_count ,
        filepath = fileparts(scopeloc.filepath{tile_index}) ;
        vecfield_path{tile_index} = filepath(length(vecfield_root)+1:end) ;  % e.g. '/2020-09-25/00/00000'
    end
    
    % Compute the baseline affine transform for each tile
    baseline_affine_transform_from_tile_index = zeros(3, 4, tile_count) ;  % nm
    for tile_index = 1 : tile_count ,
        % form an affine matrix
        linear_transform = linear_transform_from_tile_index(:,:,tile_index) ;  % nm
        nominal_tile_offset_in_mm = scopeloc.loc(tile_index,:) ;  % 1x3, mm
        nominal_tile_offset_in_nm = 1e6 * nominal_tile_offset_in_mm ;  % 1x3, nm        
        baseline_affine_transform = [linear_transform nominal_tile_offset_in_nm'] ;  % nm
        baseline_affine_transform_from_tile_index(:,:,tile_index) = baseline_affine_transform ;
    end

    % Collect some information about the point correspondences for each tile
    has_zero_z_face_matches_from_tile_index = false(tile_count,1) ;
    match_statistics = nan(tile_count, 8) ;
    for tile_index = 1 : tile_count ,
        % pix stats
        this_tile_regpts = regpts{tile_index} ;
        match_coords = this_tile_regpts.X ;  % match_count x 3, the coordinates of matched landmarks in this tile
        if isempty(match_coords) ,
            has_zero_z_face_matches_from_tile_index(tile_index) = true ;
            continue
        end
        match_coords_in_neighbor = this_tile_regpts.Y ;  % match_count x 3, the coordinates of matched landmarks in the z+1 tile, in same order as layer
        match_z_in_both_tiles = [ match_coords(:,3) match_coords_in_neighbor(:,3) ] ;  % match_count x 2
        median_match_z_in_both_tiles = round(median(match_z_in_both_tiles,1)) ; % 1x2
        min_match_z_in_both_tiles = round(min(match_z_in_both_tiles,[],1)) ; % 1x2
        max_match_z_in_both_tiles = round(max(match_z_in_both_tiles,[],1)) ; % 1x2
        tile_index_of_z_plus_1_tile = this_tile_regpts.neigs(4) ;
        match_statistics(tile_index,:) = ...
            [ tile_index  tile_index_of_z_plus_1_tile ...
              median_match_z_in_both_tiles min_match_z_in_both_tiles max_match_z_in_both_tiles ] ;
    end   
        
    % Compute the first pass as the control point targets in the rendered space
    baseline_targets_from_tile_index = zeros(cpg_ij_count*cpg_k_count, 3, tile_count) ;  % Location of each control point in the rendered space, for each tile
    for tile_index = 1 : tile_count ,
        % field curvature
        if do_apply_field_correction ,
            this_tile_curve_model = curvemodel(:,:,tile_index) ;
            field_corrected_cpg_ij1s = util.fcshift(this_tile_curve_model, order, tile_ij1s, tile_shape_ijk, cpg_ij1s) ;
                % 25 x 2, one-based ij coordinates, but non-integral
        else
            field_corrected_cpg_ij1s = cpg_ij1s ;
        end
        field_corrected_cpg_ij0s = field_corrected_cpg_ij1s - 1 ; % 25 x 2, zero-based ij coordinates, but non-integral

        % Get the affine transform
        baseline_affine_transform = baseline_affine_transform_from_tile_index(:,:,tile_index) ;
        
        % Compute the initial targets at each z plane of the control point grid
        targets_from_cpg_k_index = zeros(cpg_ij_count, 3, cpg_k_count) ;        
        for cpg_k_index = 1 : cpg_k_count ,
            field_corrected_cpg_ijk0s = zeros(cpg_ij_count, 3) ;
            field_corrected_cpg_ijk0s(:,1:2) = field_corrected_cpg_ij0s ;
            field_corrected_cpg_ijk0s(:,3) = default_cpg_k0_values(cpg_k_index) ;
            targets_at_this_cpg_k_index = add_ones_column(field_corrected_cpg_ijk0s) * baseline_affine_transform' ;
            targets_from_cpg_k_index(:, :, cpg_k_index) = targets_at_this_cpg_k_index ;
        end                   

        % Stuff it all into the baseline_targets_from_tile_index array
        %targets = reshape(targets_from_cpg_k_index, [cpg_ij_count*cpg_k_count 3]) ;
        %  ^ this gets the order wrong        
        targets = zeros(cpg_ij_count*cpg_k_count,3) ;
        offset = 0 ;
        for cpg_k_index = 1 : cpg_k_count ,
            targets(offset+1:offset+cpg_ij_count,:) = targets_from_cpg_k_index(:,:,cpg_k_index) ;
            offset = offset + cpg_ij_count ;
        end        
        baseline_targets_from_tile_index(:, :, tile_index) = targets ;
    end
    
    % The first pass of the z planes of the control points, for each tile
    baseline_cpg_k0_values_from_tile_index = repmat(default_cpg_k0_values', [1 tile_count]) ;
    
    % Sort out which we'll run with which thing
    % The tiles we run with optimpertile are "anchors".
    % The tiles we run with nomatchoptim are "floaters".
    is_floater_from_tile_index = has_k_plus_1_tile_from_tile_index & has_zero_z_face_matches_from_tile_index ;
    is_anchor_from_tile_index = ~is_floater_from_tile_index ;
    
    % Make an interpolator for each z layer, use it to shift the control point
    % targets for the tiles in the layer, and the z+1 layer
    tile_k_from_tile_index = tile_ijk_from_tile_index(:,3) ;  % column
    tile_k_from_layer_index = unique(tile_k_from_tile_index) ;  % layer of tiles, that is
    if nargin<6 || isempty(tile_k_from_run_layer_index) ,
        tile_k_from_run_layer_index = tile_k_from_layer_index(1:end-1)' ;  % a "run layer" is a layer that will actually be run
    end
    cpg_k0_values_from_tile_index = baseline_cpg_k0_values_from_tile_index ;
    targets_from_tile_index = baseline_targets_from_tile_index ;
    run_layer_count = length(tile_k_from_run_layer_index) ;
    for run_layer_index = 1 : run_layer_count ,
        % Sort out which tiles are in this layer
        tile_k = tile_k_from_run_layer_index(run_layer_index) ;
        fprintf('    Layer %d of %d, tile k/z = %d\n', run_layer_index, run_layer_count, tile_k);
        is_in_this_layer_from_tile_index = (tile_k_from_tile_index'==tile_k) ;
        tile_index_from_layer_tile_index = find(is_in_this_layer_from_tile_index);
        if isempty(tile_index_from_layer_tile_index) ,
            fprintf('No tiles found in layer with tile k/z = %d!!\n', tile_k) ;
            continue
        end
        
        % get interpolants based on paired descriptors
        [Fx_layer, Fy_layer, Fz_layer, Fx_next_layer, Fy_next_layer, Fz_next_layer, XYZ_original, XYZ_neighbor_original, outliers] =...
            util.getInterpolants(tile_index_from_layer_tile_index, ...
                                 regpts, ...
                                 baseline_affine_transform_from_tile_index, ...
                                 tile_ij1s, ...
                                 params, ...
                                 curvemodel, ...
                                 do_apply_field_correction) ;

        % Show some debugging output if called for
        if params.debug ,
            vector_field_3d_debug_script(scopeloc, scopeparams, params, XYZ_original, XYZ_neighbor_original, outliers) ;
        end
        
        % If too few matched landmarks, don't proceed with this layer
        if isempty(Fx_layer) || size(Fx_layer.Points,1) < 10 ,
            fprintf('    Layer with k/z = %d has too few matches to proceed.\n', tile_k) ;            
            continue
        end
        
        % Print the number of matches in this layer
        layer_used_match_count = size(Fx_layer.Points,1) ;
        fprintf('    Layer with k/z = %d total used matches: %d\n', tile_k, layer_used_match_count) ;
        
        % Go through all the tiles which either have no z+1 tile or for which we have
        % nonzero z-face matches, and optimize the CPG for each.  Then, compute the
        % target values for the moved z-planes of the CPG, using the interpolators to
        % try to make the matched pairs each go to the same point.
        is_anchor_from_layer_tile_index = is_anchor_from_tile_index(tile_index_from_layer_tile_index) ;
        tile_index_from_layer_anchor_index = tile_index_from_layer_tile_index(is_anchor_from_layer_tile_index) ;
        for tile_index = tile_index_from_layer_anchor_index ,
            neighbor_tile_index = tileneighbors(tile_index, 7) ;  % the z+1 tile

            if do_apply_field_correction ,
                this_tile_curve_model = curvemodel(:,:,tile_index) ;
                field_corrected_cpg_ij1s = util.fcshift(this_tile_curve_model, order, tile_ij1s, tile_shape_ijk, cpg_ij1s) ;
            else
                field_corrected_cpg_ij1s = cpg_ij1s ;
            end
            field_corrected_cpg_ij0s = field_corrected_cpg_ij1s - 1 ;

            [targets_at_cpg_k_indices_3_and_4, ...
             other_tile_targets_at_cpg_k_indices_1_and_2, ...
             new_cpg_k0_values, ...
             k_plus_1_tile_cpg_k0_values] = ...
                optimpertile(tile_index, ...
                             params, ...
                             neighbor_tile_index, ...
                             baseline_affine_transform_from_tile_index, ...
                             match_statistics, ...
                             cpg_k0_values_from_tile_index, ...
                             field_corrected_cpg_ij0s, ...
                             Fx_layer, ...
                             Fy_layer, ...
                             Fz_layer, ...
                             Fx_next_layer, ...
                             Fy_next_layer, ...
                             Fz_next_layer, ...
                             default_cpg_k0_values) ;
            cpg_k0_values_from_tile_index(:,tile_index) = new_cpg_k0_values ;
            if ~isnan(neighbor_tile_index) ,
                cpg_k0_values_from_tile_index(:,neighbor_tile_index) = k_plus_1_tile_cpg_k0_values ;
            end
            targets_from_tile_index(2*cpg_ij_count+1:end,:,tile_index) = targets_at_cpg_k_indices_3_and_4 ;
            if ~isempty(other_tile_targets_at_cpg_k_indices_1_and_2) ,
                targets_from_tile_index(1:2*cpg_ij_count,:,neighbor_tile_index) = other_tile_targets_at_cpg_k_indices_1_and_2 ;
            end
        end  % for tile_index = ...
        
        % Get the tile ijk indices of the anchors and the floaters
        is_floater_from_tile_within_layer_index = is_floater_from_tile_index(tile_index_from_layer_tile_index) ;
        tile_index_from_layer_floater_index = tile_index_from_layer_tile_index(is_floater_from_tile_within_layer_index) ;
        layer_floater_count = length(tile_index_from_layer_floater_index) ;
        tile_ijk_from_layer_anchor_index = tile_ijk_from_tile_index(tile_index_from_layer_anchor_index,:) ;        
        tile_ijk_from_layer_floater_index = tile_ijk_from_tile_index(tile_index_from_layer_floater_index,:) ;
        
        % Build an index of the closest anchor tile to each floater tile
        nearest_layer_anchor_index_from_layer_floater_index = knnsearch(tile_ijk_from_layer_anchor_index,tile_ijk_from_layer_floater_index,'K',1) ;
        
        % Build an index of all the anchor tiles in the "18-neighborhood" of each
        % floater tile.  If that set is empty for any floater tile, augment with
        % whatver the nearest anchor tile is.
        nearby_layer_anchor_indices_from_layer_floater_index = rangesearch(tile_ijk_from_layer_anchor_index,tile_ijk_from_layer_floater_index,sqrt(2)) ;
            % "nearby" here seems to mean not the 6-neighborhood nor the 26-neighborhood, but
            % the 18-neighborhood: All the neighbors that differ by 1 in at most 2
            % dimensions.  So the whole rubick's cube (3^3==27), minus self (-1), minus the
            % corners (-8).  So 27-1-8==18, the '18-neighborhood'.
        for layer_floater_index = 1 : layer_floater_count ,
            layer_anchor_index_from_nearby_tile_index = nearby_layer_anchor_indices_from_layer_floater_index{layer_floater_index} ;
            if isempty(layer_anchor_index_from_nearby_tile_index) ,
                nearby_layer_anchor_indices_from_layer_floater_index{layer_floater_index} = ...
                    nearest_layer_anchor_index_from_layer_floater_index(layer_floater_index);
            end
        end
        
        % For each floater tile, optimize its CPG and recompute targets, if needed
        for layer_floater_index = 1 : layer_floater_count ,
            % Get the tile index, and the z+1 tile index 
            tile_index = tile_index_from_layer_floater_index(layer_floater_index) ;            
            neighbor_tile_index = tileneighbors(tile_index,7) ;  % the z/k+1 tile
            
            % Get a fallback value for the CPG k0/z levels, by taking the median of the the
            % CPG k0 values for the nearby anchor tiles.
            layer_anchor_index_from_nearby_tile_index = nearby_layer_anchor_indices_from_layer_floater_index{layer_floater_index} ;            
            tile_index_from_nearby_tile_index = tile_index_from_layer_anchor_index(layer_anchor_index_from_nearby_tile_index) ;
            cpg_k0_values_from_nearby_tile_index = cpg_k0_values_from_tile_index(:,tile_index_from_nearby_tile_index) ;
            local_default_cpg_k0_values = round(median(cpg_k0_values_from_nearby_tile_index,2))' ;

            [targets_at_cpg_k_indices_3_and_4, ...
             other_tile_targets_at_cpg_k_indices_1_and_2, ...
             new_cpg_k0_values, ...
             k_plus_1_tile_cpg_k0_values] = ...
                nomatchoptim(tile_index, ...
                             params, ...
                             neighbor_tile_index, ...
                             baseline_affine_transform_from_tile_index, ...
                             match_statistics, ...
                             cpg_k0_values_from_tile_index(:,tile_index), ...
                             cpg_k0_values_from_tile_index(:,neighbor_tile_index), ...
                             field_corrected_cpg_ij0s, ...
                             Fx_layer, ...
                             Fy_layer, ...
                             Fz_layer, ...
                             Fx_next_layer, ...
                             Fy_next_layer, ...
                             Fz_next_layer, ...
                             local_default_cpg_k0_values) ;
            cpg_k0_values_from_tile_index(:,tile_index) = new_cpg_k0_values ;
            cpg_k0_values_from_tile_index(:,neighbor_tile_index) = k_plus_1_tile_cpg_k0_values ;
            targets_from_tile_index(2*cpg_ij_count+1:end,:,tile_index) = targets_at_cpg_k_indices_3_and_4 ;
            targets_from_tile_index(1:2*cpg_ij_count,:,neighbor_tile_index) = other_tile_targets_at_cpg_k_indices_1_and_2 ;
        end
    end  % for loop over layers
    
    % Compute bounding boxes for each tile
    bbox = zeros(8,3,tile_count) ;
    origin = zeros(tile_count,3) ; 
    sz = zeros(tile_count,3) ;
    for tile_index = 1 : tile_count ,
        targets = targets_from_tile_index(:,:,tile_index) ;
        [bbox(:,:,tile_index), origin(tile_index,:), sz(tile_index,:)] = ...
            util.bboxFromCorners(targets) ;
    end
    
    % Update per-tile affine transforms based on control points
    fprintf('Updating per-tile affine transforms based on control points...\n')
    final_affine_transform_from_tile_index = zeros(5, 5, tile_count) ;
    numX = size(targets_from_tile_index(:,:,1),1) ;
    for tile_index = 1 : tile_count ,
        X = targets_from_tile_index(:,:,tile_index)'/1000 ;  % /1000 converts nm to um (why? numerical stability?)
        x = cpg_i0_values ;
        y = cpg_j0_values ;
        z = cpg_k0_values_from_tile_index(:,tile_index)' ;
        [xx,yy,zz] = ndgrid(x,y,z) ;
        r = [xx(:),yy(:),zz(:)]' ;
        X_aug = [X;ones(1,numX)] ;
        r_aug = [r;ones(1,numX)] ;
        %A = (X*x')/(x*x');
        A = X_aug/r_aug ;  % 4x4, Compute a single transform that approximates, as best it can, the mapping implied by all the control-point to target pairs.
        A_full = eye(5) ;
        A_full(1:3,1:3) = A(1:3,1:3)*1000 ;  % *1000 converts um back to nm (but why convert to um in the first place?)
        A_full(1:3,5) = A(1:3,4)*1000 ;
        A_full(5,1:3) = A(4,1:3)*1000 ;  % probably not necessary, these are all zero close to zero
        final_affine_transform_from_tile_index(:,:,tile_index) = A_full ;
    end
    fprintf('Done updating per-tile affine transforms based on control points.\n')
    
    % Package everything up in a single struct for return
    vecfield = struct() ;
    vecfield.root = vecfield_root ;
    vecfield.path = vecfield_path ;
    vecfield.control = targets_from_tile_index;
    vecfield.xlim_cntrl = cpg_i0_values;
    vecfield.ylim_cntrl = cpg_j0_values;
    vecfield.zlim_cntrl = cpg_k0_values_from_tile_index;
    vecfield.afftile = baseline_affine_transform_from_tile_index;
    vecfield.tform = final_affine_transform_from_tile_index;
    vecfield.corrctrlpnttmp = field_corrected_cpg_ij0s;
    vecfield.ctrlpnttmp = cpg_ij1s-1;
    vecfield.bbox = bbox;
    vecfield.origin = origin;
    vecfield.sz = sz;
    vecfield.time = datestr(now);
    vecfield.theselayers = tile_k_from_run_layer_index;
end
