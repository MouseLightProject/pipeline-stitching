function stitch_2019_09_06()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-09-06/Tiling' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-09-06' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-09-06/stitching-output_TF4'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
