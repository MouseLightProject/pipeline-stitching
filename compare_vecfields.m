function are_close_enough = compare_vecfields(test_vecfield, reference_vecfield)
    if ~isequal(test_vecfield.xlim_cntrl, reference_vecfield.xlim_cntrl) ,
        fprintf('xlim_cntrl''s are not equal') ;
        are_close_enough = false ;
    elseif ~isequal(test_vecfield.ylim_cntrl, reference_vecfield.ylim_cntrl) ,
        fprintf('ylim_cntrl''s are not equal') ;
        are_close_enough = false ;
    elseif ~isequal(test_vecfield.zlim_cntrl, reference_vecfield.zlim_cntrl) ,
        fprintf('zlim_cntrl''s are not equal') ;
        are_close_enough = false ;
    else
        dcontrol = test_vecfield.control - reference_vecfield.control ;  % nm
        abs_diff = abs(dcontrol) ;
        max_abs_diff = max(max(abs_diff))
        median_abs_diff = median(abs_diff(:)) 
        rms_diff = sqrt( mean( (abs_diff(:)).^2 ) ) 
        are_close_enough = ( max_abs_diff < 10 ) ;
    end     
end
