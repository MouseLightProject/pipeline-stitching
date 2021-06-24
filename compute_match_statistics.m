function [total_ssd, ssd_from_tile_index, match_count_from_tile_index, rmse_from_tile_index, rmse_from_tile_ijk1] = ...
        compute_match_statistics(distance_from_match_index_from_neighbor_index_from_tile_index, ...
                                 is_in_cpg_from_match_index_from_neighbor_index_from_tile_index, ...
                                 tile_ijk1_from_tile_index, ...
                                 tile_index_from_tile_ijk1)

    % If no 2nd arg, assume all distances are valid                         
    if ~exist('is_in_cpg_from_match_index_from_neighbor_index_from_tile_index', 'var') || ...
            isempty(is_in_cpg_from_match_index_from_neighbor_index_from_tile_index) ,
        is_in_cpg_from_match_index_from_neighbor_index_from_tile_index = ...
            cellfun(@(x)(true(size(x))), ...
                    distance_from_match_index_from_neighbor_index_from_tile_index, ...
                    'UniformOutput', false) ;
    end
        
    ssd_from_neighbor_index_from_tile_index = ...
        cellfun(@sum_of_squares_after_filter, ...
                distance_from_match_index_from_neighbor_index_from_tile_index, ...
                is_in_cpg_from_match_index_from_neighbor_index_from_tile_index) ;
    match_count_from_neighbor_index_from_tile_index = cellfun(@sum, is_in_cpg_from_match_index_from_neighbor_index_from_tile_index) ;
                                           
    % sum across neighbors                                       
    ssd_from_tile_index = sum(ssd_from_neighbor_index_from_tile_index, 1) ;
    match_count_from_tile_index = sum(match_count_from_neighbor_index_from_tile_index, 1) ;

    % RMSE
    rmse_from_tile_index = sqrt(ssd_from_tile_index./match_count_from_tile_index) ;
    
    % Total SSD
    total_ssd = sum(ssd_from_tile_index) ;
    
    % make an RMSE stack
    rmse_from_tile_ijk1 = nan(size(tile_index_from_tile_ijk1)) ;
    tile_count = size(tile_ijk1_from_tile_index, 1) ;
    for tile_index = 1 : tile_count ,
        tile_ijk1 = tile_ijk1_from_tile_index(tile_index,:) ;
        rmse_from_tile_ijk1(tile_ijk1(1), tile_ijk1(2), tile_ijk1(3)) = rmse_from_tile_index(tile_index) ;
    end
    
end



function result = sum_of_squares_after_filter(distance_from_match_index, is_within_cpg_from_match_index)
    distance_from_in_cpg_match_index = distance_from_match_index(is_within_cpg_from_match_index) ;
    result = sum(distance_from_in_cpg_match_index.^2) ;    
end
