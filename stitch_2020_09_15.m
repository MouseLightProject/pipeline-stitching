sample_date = '2020-09-15'  %#ok<NOPTS>
tile_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling', sample_date)  %#ok<NOPTS>
pipeline_output_folder = sprintf('/nrs/mouselight/pipeline_output/%s', sample_date)  %#ok<NOPTS> 
stitching_output_folder_path = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/%s/stitching-output', sample_date)  %#ok<NOPTS>

stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
