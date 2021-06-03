% Specify the input/output folders
sample_date = '2020-09-15'  %#ok<NOPTS>
do_perform_field_correction = true 
tile_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling', sample_date)  %#ok<NOPTS>
pipeline_output_folder = sprintf('/nrs/mouselight/pipeline_output/%s', sample_date)  %#ok<NOPTS> 
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
memo_folder_path = fullfile(this_folder_path, sprintf('memos-%s', sample_date)) ;
stitching_output_folder_path = fullfile(memo_folder_path, 'stitching-output') ;  

% Call the function that does the real work
stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path, do_perform_field_correction) ;
