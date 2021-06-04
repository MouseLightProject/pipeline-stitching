% cd to the folder containing this file, set up the path properly
this_file_path = mfilename('fullpath') ;
this_folder_path = fileparts(this_file_path) ;
cd(this_folder_path) ;
modpath ;

% Specify the input/output folders
sample_date = '2021-05-22'  %#ok<NOPTS>
do_force_computations = false
do_perform_field_correction = true
do_show_visualizations = false
tile_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/acquisition/%s/', sample_date)  %#ok<NOPTS>
pipeline_output_folder = sprintf('/nrs/mouselight/pipeline_output/%s', sample_date)  %#ok<NOPTS> 
stitching_output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/stitching-output', sample_date)  %#ok<NOPTS>

% Call the function that does the real work
stitch(tile_folder_path, ...
       pipeline_output_folder, ...
       stitching_output_folder_path, ...
       do_force_computations, ...
       do_perform_field_correction, ...
       do_show_visualizations) ;
