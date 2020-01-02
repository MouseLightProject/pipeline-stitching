function stitch_2019_11_11_last_plane()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2019-11-12' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-11-11-last-plane' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019_11_11-last-plane/stitching-output_TF1'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
