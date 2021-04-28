function are_close_enough = compare_vecfields(test_vecfield, reference_vecfield)
    % Look at the xlims
    test_xlim_cntrl = test_vecfield.xlim_cntrl ;
    reference_xlim_cntrl = reference_vecfield.xlim_cntrl ;
    diff_xlim_cntrl = test_xlim_cntrl - reference_xlim_cntrl
    xlim_close_enough = all(diff_xlim_cntrl==0)
    
    % Look at the ylims
    test_ylim_cntrl = test_vecfield.ylim_cntrl ;
    reference_ylim_cntrl = reference_vecfield.ylim_cntrl ;
    diff_ylim_cntrl = test_ylim_cntrl - reference_ylim_cntrl
    ylim_close_enough = all(diff_ylim_cntrl==0)
    
    % Look at the zlims
    test_zlim_cntrl = test_vecfield.zlim_cntrl ;
    reference_zlim_cntrl = reference_vecfield.zlim_cntrl ;
    diff_zlim_cntrl = test_zlim_cntrl - reference_zlim_cntrl ;
    do_zlims_differ_from_tile_index = any(diff_zlim_cntrl>0) ;
    tile_index_from_different_index = find(do_zlims_differ_from_tile_index)
    differing_test_zlim_cntrl = test_zlim_cntrl(:,tile_index_from_different_index)
    differing_reference_zlim_cntrl = reference_zlim_cntrl(:,tile_index_from_different_index)
    zlim_close_enough = length(tile_index_from_different_index)<=5 && max(max(abs(diff_zlim_cntrl)))<=2 

    % Look at the control points for tiles with the same zlims
    %are_zlims_same_from_tile_index = ~do_zlims_differ_from_tile_index ;    
    raw_dcontrol = test_vecfield.control - reference_vecfield.control ;  % nm
    dcontrol = raw_dcontrol ;
    dcontrol(:,:,do_zlims_differ_from_tile_index) = 0 ;  % differences in these ones are to be expected, and we don't really care about them
    abs_diff = abs(dcontrol) ;
    [max_abs_diff, tile_index_of_max_abs_diff] = max(max(max(abs_diff)))
    
    %test_vecfield.control(:,:,tile_index_of_max_abs_diff)
    %reference_vecfield.control(:,:,tile_index_of_max_abs_diff)    
    
    median_abs_diff = median(abs_diff(:)) 
    rms_diff = sqrt( mean( (abs_diff(:)).^2 ) ) 
    are_targets_close_enough = ( max_abs_diff < 50 ) ;
    
    are_close_enough = xlim_close_enough && ylim_close_enough && zlim_close_enough && are_targets_close_enough ;
end
