% Specify the input/output folders
sample_date = '2020-09-15'  %#ok<NOPTS>
do_force_computations = false 
do_perform_field_correction = true 
do_show_visualizations = false 
tile_folder_path = sprintf('/nearline/mouselight/data/RAW_archive/%s/Tiling', sample_date)  %#ok<NOPTS>
pipeline_output_folder = sprintf('/nrs/mouselight/pipeline_output/%s', sample_date)  %#ok<NOPTS> 
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
memo_folder_path = fullfile(this_folder_path, sprintf('memos-%s', sample_date)) ;
stitching_output_folder_path = fullfile(memo_folder_path, 'stitching-output') ;  

% Call the function that does the real work
tic_id = tic() ;
stitch(tile_folder_path, ...
       pipeline_output_folder, ...
       stitching_output_folder_path, ...
       do_force_computations, ...
       do_perform_field_correction, ...
       do_show_visualizations) ;

vecfield_file_path = fullfile(stitching_output_folder_path, 'vecfield3D.mat')
mat = load(vecfield_file_path, 'vecfield3D', 'params') ;
vecfield = mat.vecfield3D ;

parent_folder_path = fileparts(this_folder_path) ;
reference_vecfield_file_path = fullfile(parent_folder_path, 'pipeline-stitching-fb8d827/2020-09-15-stitching-output-verify/vecfield3D.mat')
reference_mat = load(reference_vecfield_file_path, 'vecfield3D', 'params') ;
reference_vecfield = reference_mat.vecfield3D ;

are_close_enough = compare_vecfields(vecfield, reference_vecfield)
elapsed_time = toc(tic_id) ;
fprintf('Elapsed time was %g seconds.\n', elapsed_time) ;
