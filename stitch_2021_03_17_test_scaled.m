% cd to the folder containing this file, set up the path properly
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
cd(this_folder_path) ;
modpath ;

% Specify the input/output folders
tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/test-data/test-scaled/2021-03-17'
pipeline_output_folder = '/nrs/mouselight/pipeline_output/test-scaled/2021-03-17'
stitching_output_folder_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/test-scaled/2021-03-17/stitching-output'

% Call the function that does the real work
stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
