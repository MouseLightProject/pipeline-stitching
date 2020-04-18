tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-10-04/Tiling'  %#ok<NOPTS>
pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-10-04' %#ok<NOPTS>
stitching_output_folder_path = '/groups/mousebrainmicro/mousebrainmicro/cluster/Reconstructions/2019-10-04/stitching-output-alt-test'  %#ok<NOPTS>

stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
