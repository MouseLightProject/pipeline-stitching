function vecfield = vectorField3D(params, scopeloc, regpts, scopeparams, curvemodel, tile_k_from_run_layer_index)
    fprintf('Calculating vector fields\n') ;
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
    do_apply_FC = 1 ;  % whether or not to apply the field correction, I think

    % Build the neighbor index
    tileneighbors = buildNeighbor(scopeloc.gridix(:,1:3)) ;  % tile_count x 7, each row in [self -x -y +x +y -z +z] format
    
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
    
    % Compute the affine transform for each tile
    affine_transform_from_tile_index = zeros(3, 4, tile_count) ;  % nm
    for tile_index = 1 : tile_count ,
        % form an affine matrix
        linear_transform = linear_transform_from_tile_index(:,:,tile_index) ;  % nm
        nominal_tile_offset_in_mm = scopeloc.loc(tile_index,:) ;  % 1x3, mm
        nominal_tile_offset_in_nm = 1e6 * nominal_tile_offset_in_mm ;  % 1x3, nm        
        affine_transform = [linear_transform nominal_tile_offset_in_nm'] ;  % nm
        affine_transform_from_tile_index(:,:,tile_index) = affine_transform ;
    end

    % Compute the first pass as the control point targets in the rendered space
    targets_from_tile_index = zeros(cpg_ij_count*cpg_k_count, 3, tile_count) ;  % Location of each control point in the rendered space, for each tile
    for tile_index = 1 : tile_count ,
        % field curvature
        this_tile_curve_model = curvemodel(:,:,tile_index) ;
        field_corrected_cpg_ij1s = util.fcshift(this_tile_curve_model, order, tile_ij1s, tile_shape_ijk, cpg_ij1s) ;  
          % 25 x 2, one-based ij coordinates, but non-integral
        field_corrected_cpg_ij0s = field_corrected_cpg_ij1s - 1 ; % 25 x 2, zero-based ij coordinates, but non-integral

        % Get the affine transform
        affine_transform = affine_transform_from_tile_index(:,:,tile_index) ;
        
        % Compute the initial targets at each z plane of the control point grid
        targets_from_cpg_k_index = zeros(cpg_ij_count, 3, cpg_k_count) ;        
        for cpg_k_index = 1 : cpg_k_count ,
            field_corrected_cpg_ijk0s = zeros(cpg_ij_count, 3) ;
            field_corrected_cpg_ijk0s(:,1:2) = field_corrected_cpg_ij0s ;
            field_corrected_cpg_ijk0s(:,3) = default_cpg_k0_values(cpg_k_index) ;
            targets_at_this_cpg_k_index = add_ones_column(field_corrected_cpg_ijk0s) * affine_transform' ;
            targets_from_cpg_k_index(:, :, cpg_k_index) = targets_at_this_cpg_k_index ;
        end                   

        % Stuff it all into the targets_from_tile_index array
        %targets = reshape(targets_from_cpg_k_index, [cpg_ij_count*cpg_k_count 3]) ;
        %  ^ this gets the order wrong        
        targets = zeros(cpg_ij_count*cpg_k_count,3) ;
        offset = 0 ;
        for cpg_k_index = 1 : cpg_k_count ,
            targets(offset+1:offset+cpg_ij_count,:) = targets_from_cpg_k_index(:,:,cpg_k_index) ;
            offset = offset + cpg_ij_count ;
        end        
        targets_from_tile_index(:, :, tile_index) = targets ;
    end

    % The first pass of the z planes of the control points, for each tile
    cpg_k0_values_from_tile_index = repmat(default_cpg_k0_values', [1 tile_count]) ;
    
    % Collect some information about the point correspondences for each tile
    match_statistics = nan(tile_count, 8) ;
    for tile_index = 1 : tile_count ,
        % pix stats
        this_tile_regpts = regpts{tile_index} ;
        match_coords = this_tile_regpts.X ;  % match_count x 3, the coordinates of matched landmarks in this tile
        match_coords_in_neighbor = this_tile_regpts.Y ;  % match_count x 3, the coordinates of matched landmarks in the z+1 tile, in same order as layer
        if isempty(match_coords) ,
            continue
        end
        match_z_in_both_tiles = [ match_coords(:,3) match_coords_in_neighbor(:,3) ] ;  % match_count x 2
        median_match_z_in_both_tiles = round(median(match_z_in_both_tiles,1)) ; % 1x2
        min_match_z_in_both_tiles = round(min(match_z_in_both_tiles,[],1)) ; % 1x2
        max_match_z_in_both_tiles = round(max(match_z_in_both_tiles,[],1)) ; % 1x2
        tile_index_of_z_plus_1_tile = this_tile_regpts.neigs(4) ;
        match_statistics(tile_index,:) = ...
            [ tile_index  tile_index_of_z_plus_1_tile ...
              median_match_z_in_both_tiles min_match_z_in_both_tiles max_match_z_in_both_tiles ] ;
    end   
    
    tile_k_from_tile_index = scopeloc.gridix(:,3) ;  % column
    tile_k_from_layer_index = unique(tile_k_from_tile_index) ;  % layer of tiles, that is
    if nargin<6 || isempty(tile_k_from_run_layer_index) ,
        tile_k_from_run_layer_index = tile_k_from_layer_index(1:end-1)' ;  % a "run layer" is a layer that will actually be run
    end
    do_run_nomatchoptim_from_tile_index = NaN(1,tile_count) ;
    run_layer_count = length(tile_k_from_run_layer_index) ;
    for run_layer_index = 1 : run_layer_count ,
        % Sort out which tiles are in this layer
        tile_k = tile_k_from_run_layer_index(run_layer_index) ;
        fprintf('    Layer %d of %d, tile k/z = %d\n', run_layer_index, run_layer_count, tile_k);
        is_in_this_layer_from_tile_index = (tile_k_from_tile_index'==tile_k) ;
        tile_index_from_tile_within_layer_index = find(is_in_this_layer_from_tile_index);
        if isempty(tile_index_from_tile_within_layer_index) ,
            fprintf('No tiles found in layer with tile k/z = %d!!\n', tile_k) ;
            continue
        end
        
        % get interpolants based on paired descriptors
        [Fx_layer, Fy_layer, Fz_layer, Fx_next_layer, Fy_next_layer, Fz_next_layer, XYZ_original, XYZ_neighbor_original, outliers] =...
            util.getInterpolants(tile_index_from_tile_within_layer_index, regpts, affine_transform_from_tile_index, tile_ij1s, params, curvemodel, do_apply_FC) ;

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
        
        for tile_index = tile_index_from_tile_within_layer_index ,  % layer t
            neighbor_tile_index = tileneighbors(tile_index, 7) ;  % the z+1 tile

            this_tile_curve_model = curvemodel(:,:,tile_index) ;
            field_corrected_cpg_ij1s = util.fcshift(this_tile_curve_model, order, tile_ij1s, tile_shape_ijk, cpg_ij1s) ;
            field_corrected_cpg_ij0s = field_corrected_cpg_ij1s - 1 ;

%             if tile_index == 2137 || neighbor_tile_index == 2137 ,
%                keyboard
%             end             
            
            [do_run_nomatchoptim, ...
             control_t_bot12, ...
             control_tp1_top12, ...
             cpg_k0_values_from_tile_index] = ...
                optimpertile(tile_index, ...
                             params, ...
                             tileneighbors, ...
                             affine_transform_from_tile_index, ...
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
            do_run_nomatchoptim_from_tile_index(tile_index) = do_run_nomatchoptim ;
            if do_run_nomatchoptim ,
                % do nothing, will call nomatchoptim() for this tile below
            elseif isempty(control_tp1_top12) ,  % no below adjacent tile
                targets_from_tile_index(2*cpg_ij_count+1:end,:,tile_index) = control_t_bot12 ;
            else
                targets_from_tile_index(2*cpg_ij_count+1:end,:,tile_index) = control_t_bot12 ;
                targets_from_tile_index(1:2*cpg_ij_count,:,neighbor_tile_index) = control_tp1_top12 ;
            end

        end
        
        anchorinds = tile_index_from_tile_within_layer_index(~do_run_nomatchoptim_from_tile_index(tile_index_from_tile_within_layer_index)) ;
        anchors = scopeloc.gridix(anchorinds,:) ;        
        queryinds = tile_index_from_tile_within_layer_index(do_run_nomatchoptim_from_tile_index(tile_index_from_tile_within_layer_index)>0) ;
        queries = scopeloc.gridix(queryinds,:) ;
        IDX_nn = knnsearch(anchors,queries,'K',1) ;
        IDX = rangesearch(anchors,queries,sqrt(2)) ;
        
        for ine = 1 : length(queryinds) ,
            if isempty(IDX{ine})
                IDX{ine} = IDX_nn(ine);
            end
        end
        
        for tile_index = tile_index_from_tile_within_layer_index ,  % layer t
            neighbor_tile_index = tileneighbors(tile_index,7) ;
            
%             if tile_index == 2137 || neighbor_tile_index == 2137 ,
%                keyboard
%             end            

            if do_run_nomatchoptim_from_tile_index(tile_index) ,
                is_query =(tile_index==queryinds) ;
                local_default_cpg_k0_values = round(median(cpg_k0_values_from_tile_index(:,anchorinds(IDX{is_query})),2))' ;
                
                [do_run_nomatchoptim_from_tile_index(tile_index), ...
                 control_t_bot12, ...
                 control_tp1_top12, ...
                 cpg_k0_values_from_tile_index] = ...
                    nomatchoptim(tile_index, ...
                                 params, ...
                                 tileneighbors, ...
                                 affine_transform_from_tile_index, ...
                                 match_statistics, ...
                                 cpg_k0_values_from_tile_index, ...
                                 field_corrected_cpg_ij0s, ...
                                 Fx_layer, ...
                                 Fy_layer, ...
                                 Fz_layer, ...
                                 Fx_next_layer, ...
                                 Fy_next_layer, ...
                                 Fz_next_layer, ...
                                 local_default_cpg_k0_values) ;

                targets_from_tile_index(2*cpg_ij_count+1:end,:,tile_index) = control_t_bot12 ;
                targets_from_tile_index(1:2*cpg_ij_count,:,neighbor_tile_index) = control_tp1_top12 ;
            end
        end
    end
    
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
    affine_transform = zeros(5,5,tile_count);
    numX = size(targets_from_tile_index(:,:,1),1);
    for tile_index = 1 : tile_count ,
        X = targets_from_tile_index(:,:,tile_index)'/1000;
        x = cpg_i0_values;
        y = cpg_j0_values;
        z = cpg_k0_values_from_tile_index(:,tile_index)';
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
        affine_transform(:,:,tile_index) = Aest;
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
    vecfield.afftile = affine_transform_from_tile_index;
    vecfield.tform = affine_transform;
    vecfield.corrctrlpnttmp = field_corrected_cpg_ij0s;
    vecfield.ctrlpnttmp = cpg_ij1s-1;
    vecfield.bbox = bbox;
    vecfield.origin = origin;
    vecfield.sz = sz;
    vecfield.time = datestr(now);
    vecfield.theselayers = tile_k_from_run_layer_index;
end
