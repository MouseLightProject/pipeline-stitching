function stitch_2019_10_04()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-10-04/Tiling' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-10-04' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-10-04/stitching-output_TF1'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
