% Specify the input/output folders
sample_date = '2021-03-17'  %#ok<NOPTS>
tile_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling', sample_date)  %#ok<NOPTS>
pipeline_output_folder = sprintf('/nrs/mouselight/pipeline_output/%s', sample_date)  %#ok<NOPTS> 
this_folder_path = fileparts(mfilename('fullpath')) ;
stitching_output_folder_path = fullfile(this_folder_path, sprintf('%s-stitching-output-alt', sample_date))  %#ok<NOPTS>

% Call the function that does the real work
stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
