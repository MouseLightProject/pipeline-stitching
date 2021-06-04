% cd to the folder containing this file, set up the path properly
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
cd(this_folder_path) ;
modpath ;

% Specify the input/output folders
tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/test-data/test-unscaled/2021-03-17'
pipeline_output_folder = '/nrs/mouselight/pipeline_output/test-unscaled/2021-03-17-new-classifier'
stitching_output_folder_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/test-unscaled/2021-03-17-new-classifier/stitching-output'

% Call the function that does the real work
stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
