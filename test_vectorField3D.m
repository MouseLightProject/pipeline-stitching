tic_id = tic() ;
this_folder_path = fileparts(mfilename('fullpath')) ;
data_file_path = fullfile(this_folder_path, 'vectorfield3D-test-inputs-from-2020-09-15.mat') ;
load(data_file_path, 'params', 'scopeloc', 'regpts', 'scopeparams', 'curvemodel', 'tile_k_from_run_layer_index') ;
vecfield = vectorField3D(params, scopeloc, regpts, scopeparams, curvemodel, tile_k_from_run_layer_index) ;

parent_folder_path = fileparts(this_folder_path) ;
reference_vecfield_file_path = fullfile(parent_folder_path, 'pipeline-stitching-fb8d827/2020-09-15-stitching-output-verify/vecfield3D.mat')
reference_mat = load(reference_vecfield_file_path, 'vecfield3D', 'params') ;
reference_vecfield = reference_mat.vecfield3D ;

are_close_enough = compare_vecfields(vecfield, reference_vecfield)
elapsed_time = toc(tic_id) ;
fprintf('Elapsed time was %g seconds.\n', elapsed_time) ;
