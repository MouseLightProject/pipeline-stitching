% Specify the input/output folders
sample_tag = '2021-03-17-test-unscaled'  %#ok<NOPTS>
analysis_tag = 'adam-classifier' 
do_force_computations = true
do_perform_field_correction = true
do_show_visualizations = true
this_folder_path = fileparts(mfilename('fullpath')) ;
tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/test-data/test-unscaled/2021-03-17'
pipeline_output_folder = '/nrs/mouselight/pipeline_output/test-unscaled/2021-03-17-adam-classifier'
sample_memo_folder_path = fullfile(this_folder_path, sprintf('memos-%s', sample_tag)) ;
analysis_memo_folder_path = fullfile(sample_memo_folder_path, analysis_tag) ;
stitching_output_folder_path = fullfile(analysis_memo_folder_path, 'stitching-output') ;

% Call the function that does the real work
stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path, do_force_computations, do_perform_field_correction, do_show_visualizations) ;
