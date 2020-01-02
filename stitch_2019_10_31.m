function stitch_2019_10_31()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-10-31/Tiling' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-10-31' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-10-31/stitching-output_TF1'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
