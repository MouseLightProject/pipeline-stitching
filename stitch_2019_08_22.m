function stitch_2019_08_22()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/2019-08-22/Tiling' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-08-22' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-08-22/stitching-output_TF20190914Ch1'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
