function xyz_from_query_index = ...
        make_and_apply_barycentric_interpolant(cpg_i0_values , ...
                                               cpg_j0_values , ...
                                               cpg_k0_values , ...
                                               xyz_from_grid_index, ...
                                               ijk0_from_query_index)
    % "CPG" stands for "control point grid".  It's the grid of input points that we
    % use to constuct the interpolator.
    % cpg_i0_values: row vector of x/i values in the CPG
    % cpg_j0_values: row vector of y/j values in the CPG
    % cpg_k0_values: row vector of z/k values in the CPG
    % xyz_from_grid_index: For each point on the CPG, the point in 3-space that
    %   it should be mapped to, grid_count x 3, where grid_count is the product of
    %   the lengths of cpg_i0_values, cpg_j0_values, and cpg_k0_values.  These must
    %   be in "z-major" order.  I.e. x coord changes most quickly, y next, then z.
    % ijk0_from_query_index, the x/i, y/j, and z/k coords of each query point to be fed
    %   into the interpolator. query_count x 3
    %
    % Output:
    %   xyz_from_query_index, the output xyz point for each query.  query_count x 3
        
    % Make the interpolator
    bi = ...
        barycentric_interpolant(cpg_i0_values , ...
                                cpg_j0_values , ...
                                cpg_k0_values , ...
                                xyz_from_grid_index) ;
    
    % Use the interpolator                        
    xyz_from_query_index = bi(ijk0_from_query_index) ;                            
end
