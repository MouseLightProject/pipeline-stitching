function stitch_2019_11_11_last_plane_alt()    
    tile_folder_path = '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2019-11-12' %#ok<NOPRT>
    pipeline_output_folder =  '/nrs/mouselight/pipeline_output/2019-11-11-last-plane' %#ok<NOPRT>
    stitching_output_folder_path = '/nrs/mouselight/cluster/classifierOutputs/2019-11-11-last-plane/stitching-output-alt-1'  %#ok<NOPRT>
    stitch(tile_folder_path, pipeline_output_folder, stitching_output_folder_path) ;
end
git 