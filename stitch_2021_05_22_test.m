% Specify the input/output folders
sample_tag = '2021-05-22-test'  %#ok<NOPTS>
analysis_tag = 'production-classifier' 
do_force_computations = false
do_perform_field_correction = true
do_run_in_debug_mode = false
do_show_visualizations = true
this_folder_path = fileparts(mfilename('fullpath')) ;
tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/test-data/2021-05-22-test'
pipeline_output_folder = '/nrs/mouselight/pipeline_output/2021-05-22'
sample_memo_folder_path = fullfile(this_folder_path, 'memos', sample_tag) ;
analysis_memo_folder_path = fullfile(sample_memo_folder_path, analysis_tag) ;
stitching_output_folder_path = fullfile(analysis_memo_folder_path, 'stitching-output') ;

% Call the function that does the real work
stitch_options = struct('do_force_computations', do_force_computations, ...
                        'do_perform_field_correction', do_perform_field_correction, ...
                        'do_run_in_debug_mode', do_run_in_debug_mode, ...
                        'do_show_visualizations', do_show_visualizations) ;
stitch(tile_folder_path, ...
       pipeline_output_folder, ...
       stitching_output_folder_path, ...
       stitch_options) ;
   