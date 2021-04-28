function are_close_enough = compare_vecfields_on_disk(test_vecfield_file_name, reference_vecfield_file_name)
    reference_mat = load(reference_vecfield_file_name, 'vecfield3D', 'params') ;
    reference_vecfield = reference_mat.vecfield3D ;
    test_mat = load(test_vecfield_file_name, 'vecfield3D', 'params') ;
    test_vecfield = test_mat.vecfield3D ;    
    are_close_enough = compare_vecfields(test_vecfield, reference_vecfield) ;
end
