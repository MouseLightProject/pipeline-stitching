function [curvemodel, scopeparams] = estimate_field_curvature_for_all_tiles(scopeloc, descriptorfolder, desc_ch, params)
    if params.applyFC ,
        checkthese = [1 4 5 7];  % self - right - bottom - below
        neighbors = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format    
        % accumulator for steps
        model = @(p,y) p(3) - p(2).*((y-p(1)).^2) ;  % FC model

        % get tile descriptors
        descriptors = getDescriptorsPerFolder(descriptorfolder, scopeloc, desc_ch) ;
        fprintf('Loaded descriptors\n')

        % curvature estimation using descriptor match
        [paireddescriptor_from_xy_match, curvemodel_from_xy_match] = match.xymatch(descriptors, neighbors(:,checkthese), scopeloc, params, model) ;
        fprintf('X&Y descriptor match done\n')

        % interpolate tiles with missing parameters from adjecent tiles
        [paireddescriptor_from_COE, curvemodel_from_COE, unreliable] = ...
            match.curvatureOutlierElimination(paireddescriptor_from_xy_match, curvemodel_from_xy_match, scopeloc, params, model) ;
        fprintf('Outlier elimination done\n')

        % tile base affine
        if params.singleTile ,
            [scopeparams, curvemodel] = ...
                homographyPerTile6Neighbor(params, neighbors, scopeloc, paireddescriptor_from_COE, curvemodel_from_COE) ;
            fprintf('Per-tile affine estimation done\n')
        else
            % joint affine estimation
            scopeparams_from_EJA = ...
                match.estimatejointaffine(paireddescriptor_from_COE, neighbors, scopeloc, params, curvemodel_from_COE, 0) ;
            [scopeparams, curvemodel] = ...
                match.affineOutlierElimination( scopeloc, scopeparams_from_EJA, paireddescriptor_from_COE, curvemodel_from_COE, unreliable) ;
            fprintf('Joint affine estimation\n')
        end
    else
        % No field correction, so just use fallback values for all these things
        tile_count = length(scopeloc.relativepaths) ;
        single_tile_null_curvemodel = ...
            [ 0 0 1 ; ...
              0 0 1 ; ...
              nan nan nan ] ;
        curvemodel = repmat(single_tile_null_curvemodel, [1 1 tile_count]) ;        
        
        spacing = 1000 * params.imsize_um ./ (params.imagesize-1) ;
        default_linear_transform = diag([+1 +1 -1] .* spacing) ;
        
        scopeparams = struct_with_shape_and_fields([1 tile_count], {'imsize_um', 'dims', 'affinegl', 'affineglFC'}) ;
        for tile_index = 1 : tile_count ,       
            scopeparams(tile_index).imsize_um = params.imsize_um ;
            scopeparams(tile_index).dims = params.imagesize ;
            scopeparams(tile_index).affinegl = default_linear_transform ;
            scopeparams(tile_index).affineglFC = default_linear_transform ;            
        end
    end
end
