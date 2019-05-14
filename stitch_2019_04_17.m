function stitch_2019_04_17()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-04-17/Tiling' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-04-17' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-04-17/stitching-output'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
